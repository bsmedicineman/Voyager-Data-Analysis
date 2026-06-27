#!/usr/bin/env python3
"""
download_em_data.py
===================
Concurrently download EM (magnetometer / E-field / search-coil) CDF data for the
inner-heliosphere & Mercury missions, sort it on disk by mission/dataset, and
keep a crash-safe manifest so an interrupted run resumes exactly where it left
off.

  * discovery .... NASA CDAWeb web services (cdasws) -> exact source file URLs
  * download ..... parallel threads, chunked streaming, HTTP-Range resume
  * recovery ..... SQLite manifest + atomic .part files; re-run = resume
  * memory ....... streams to disk in small chunks; never buffers a whole file

Layout on disk:
    <output-dir>/
        download_manifest.db
        <MISSION>/<DATASET_ID>/<original_filename>.cdf

Usage
-----
  # See exactly which dataset IDs exist right now for a mission:
  python download_em_data.py list-datasets --missions psp solo

  # Download a date range for several missions, 8 parallel workers:
  python download_em_data.py fetch --missions psp solo messenger helios \\
        --start 2022-01-01 --end 2022-02-01 \\
        --output-dir /data/spacefoam --workers 8

  # Just one dataset; preview without downloading:
  python download_em_data.py fetch --missions psp \\
        --datasets PSP_FLD_L2_MAG_RTN_1MIN \\
        --start 2021-01-01 --end 2021-12-31 -o /data/sf --dry-run

  # Download from a plain list of URLs (no cdasws needed):
  python download_em_data.py fetch --url-list urls.txt -o /data/sf

  python download_em_data.py selftest      # crash-recovery logic, no network

Note: BepiColombo data lives in ESA's PSA (different protocol) and is not on
CDAWeb; those specs are flagged and skipped here with guidance. Everything else
(PSP, Solar Orbiter, MESSENGER, Helios) downloads automatically.
"""

from __future__ import annotations

import argparse
import os
import queue
import signal
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from typing import Optional

import em_foam_common as C

LOG = C.get_logger("download")
STOP = threading.Event()        # set on Ctrl-C / SIGTERM for graceful shutdown


# --------------------------------------------------------------------------- #
#  Discovery via CDAWeb
# --------------------------------------------------------------------------- #
def discover_cdaweb(specs, start: str, end: str) -> "list[tuple[str, str, Optional[int]]]":
    """Return (url, dataset_id, expected_bytes|None) for every source file."""
    try:
        from cdasws import CdasWs
    except Exception as e:
        LOG.error("cdasws not installed (%s). Use --url-list, or "
                  "`pip install cdasws`.", e)
        return []
    cdas = CdasWs()
    found: list[tuple[str, str, object]] = []
    for spec in specs:
        if spec.source != "cdaweb":
            LOG.warning("skip %s: source=%s (not automated here) -- %s",
                        spec.dataset_id, spec.source, spec.note)
            continue
        try:
            status, files = cdas.get_original_files(spec.dataset_id, start, end)
        except Exception as e:
            LOG.error("discovery failed for %s: %s", spec.dataset_id, e)
            continue
        n0 = len(found)
        for fi in (files or []):
            url = _first(fi, "Name", "name", "URL", "url")
            if not url:
                continue
            size = _first(fi, "Length", "length", "FileSize", "size")
            try:
                size = int(size) if size is not None else None
            except (TypeError, ValueError):
                size = None
            found.append((url, spec.dataset_id, size))
        LOG.info("  %-34s %4d files  %s..%s",
                 spec.dataset_id, len(found) - n0, start, end)
    return found


def _first(d, *keys):
    for k in keys:
        if isinstance(d, dict) and k in d and d[k]:
            return d[k]
        if hasattr(d, k) and getattr(d, k):
            return getattr(d, k)
    return None





def mission_of(dataset_id: str) -> str:
    for m, specs in C.MISSION_REGISTRY.items():
        if any(s.dataset_id == dataset_id for s in specs):
            return m.upper()
    return "MISC"


def local_path_for(out_dir: str, dataset_id: str, url: str) -> str:
    fname = url.split("?")[0].rstrip("/").split("/")[-1] or "file.cdf"
    return os.path.join(out_dir, mission_of(dataset_id), dataset_id, fname)


# --------------------------------------------------------------------------- #
#  Download one file (streamed, resumable)
# --------------------------------------------------------------------------- #
def download_one(session, url: str, path: str, expected, args) -> Optional[int]:
    """Stream url -> path with Range resume. Returns bytes, or None if stopped."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    part = path + ".part"
    resume_from = os.path.getsize(part) if os.path.exists(part) else 0
    headers = {"Range": f"bytes={resume_from}-"} if resume_from else {}

    with session.get(url, headers=headers, stream=True, timeout=args.timeout) as r:
        if r.status_code == 416:                      # already complete
            if os.path.exists(part):
                os.replace(part, path)
            return os.path.getsize(path)
        mode = "ab"
        if resume_from and r.status_code == 200:      # server ignored Range
            resume_from, mode = 0, "wb"
        r.raise_for_status()

        total = expected
        cl = r.headers.get("Content-Length")
        if cl:
            total = (resume_from + int(cl)) if r.status_code == 206 else int(cl)

        got = resume_from
        with open(part, mode) as fh:
            for chunk in r.iter_content(chunk_size=args.chunk_bytes):
                if STOP.is_set():
                    fh.flush()
                    return None                        # keep .part for resume
                if chunk:
                    fh.write(chunk)
                    got += len(chunk)

    if total and got < total:                          # truncated -> keep .part
        raise IOError(f"short read {got}/{total} bytes")
    os.replace(part, path)
    return got


def worker(name, work_q, manifest, session, args):
    while not STOP.is_set():
        try:
            url, dataset, path, expected = work_q.get_nowait()
        except queue.Empty:
            return
        try:
            manifest.mark(url, "downloading", bump_attempt=True)
            got = download_one(session, url, path, expected, args)
            if got is None:                            # graceful stop mid-file
                manifest.mark(url, "pending")
            else:
                manifest.mark(url, "done", got=got)
        except Exception as e:
            manifest.mark(url, "failed", error=str(e))
            LOG.debug("fail %s: %s", url, e)
        finally:
            work_q.task_done()


# --------------------------------------------------------------------------- #
#  fetch driver
# --------------------------------------------------------------------------- #
def load_url_list(path: str) -> list[tuple[str, str, object]]:
    rows = []
    with open(path) as fh:
        for ln in fh:
            u = ln.strip()
            if u and not u.startswith("#"):
                rows.append((u, "URLLIST", None))
    return rows


def fetch(args) -> int:
    os.makedirs(args.output_dir, exist_ok=True)

    if args.url_list:
        discovered = load_url_list(args.url_list)
        LOG.info("loaded %d URLs from %s", len(discovered), args.url_list)
    else:
        if not (args.missions and args.start and args.end):
            LOG.error("need --missions, --start, --end (or use --url-list)")
            return 2
        specs = C.resolve_specs(
            args.missions, set(args.datasets) if args.datasets else None)
        LOG.info("discovering source files via CDAWeb...")
        discovered = discover_cdaweb(specs, args.start, args.end)

    if not discovered:
        LOG.error("nothing to download")
        return 1

    # Resolve local paths and load the manifest.
    rows = []
    for url, dataset, size in discovered:
        path = (local_path_for(args.output_dir, dataset, url)
                if dataset not in ("URLLIST",)
                else os.path.join(args.output_dir, "URLLIST",
                                  url.split("/")[-1] or "file.cdf"))
        rows.append((url, dataset, path, size))

    if args.dry_run:
        known = sum(s for *_, s in rows if isinstance(s, int))
        LOG.info("DRY RUN: %d files, ~%s known size", len(rows),
                 C.human_bytes(known) if known else "size unknown")
        for url, ds, path, size in rows[:20]:
            LOG.info("  %s  ->  %s  (%s)", ds, os.path.relpath(path, args.output_dir),
                     C.human_bytes(size) if isinstance(size, int) else "?")
        if len(rows) > 20:
            LOG.info("  ... and %d more", len(rows) - 20)
        return 0

    manifest = C.ManifestDB(os.path.join(args.output_dir, "download_manifest.db"))
    added = manifest.upsert_pending([(u, d, p, s if isinstance(s, int) else None)
                                     for u, d, p, s in rows])
    recovered = manifest.reset_missing()
    LOG.info("manifest: +%d new, %d recovered (missing files), state=%s",
             added, recovered, manifest.stats())

    todo = manifest.get_work()
    if not todo:
        LOG.info("everything already downloaded. nothing to do.")
        manifest.close()
        return 0

    work_q: queue.Queue = queue.Queue()
    for item in todo:
        work_q.put(item)

    import requests
    session = requests.Session()
    session.headers.update({"User-Agent": "spacefoam-em/1.0"})

    LOG.info("downloading %d files with %d workers (Ctrl-C to pause+resume later)",
             len(todo), args.workers)
    t0 = time.time()

    monitor = threading.Thread(target=_progress, args=(manifest, len(todo)),
                               daemon=True)
    monitor.start()

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for i in range(args.workers):
            pool.submit(worker, f"w{i}", work_q, manifest, session, args)
        try:
            while any(t for t in [work_q.unfinished_tasks]) and not STOP.is_set():
                time.sleep(0.3)
            work_q.join()
        except KeyboardInterrupt:
            STOP.set()
            LOG.warning("stop requested -- finishing current chunks, run again to resume")

    st = manifest.stats()
    LOG.info("finished in %.1fs  state=%s", time.time() - t0, st)
    if st.get("failed"):
        LOG.info("%d files failed; just re-run the same command to retry them.",
                 st["failed"])
    manifest.close()
    return 0


def _progress(manifest: C.ManifestDB, total: int):
    while not STOP.is_set():
        time.sleep(5)
        st = manifest.stats()
        done = st.get("done", 0)
        LOG.info("progress: %d/%d done  (%d failed, %d in-flight)",
                 done, total, st.get("failed", 0), st.get("downloading", 0))
        if done + st.get("failed", 0) >= total:
            return


# --------------------------------------------------------------------------- #
#  list-datasets (confirm current CDAWeb IDs)
# --------------------------------------------------------------------------- #
def list_datasets(args) -> int:
    try:
        from cdasws import CdasWs
    except Exception as e:
        LOG.error("cdasws needed for live listing: %s", e)
        return 2
    cdas = CdasWs()
    obs_map = {"psp": "Parker Solar Probe", "solo": "Solar Orbiter",
               "messenger": "MESSENGER", "helios": "Helios"}
    for m in args.missions:
        key = m.lower()
        LOG.info("=== %s : registry defaults ===", key)
        for s in C.MISSION_REGISTRY.get(key, []):
            print(f"    {s.dataset_id:34} [{s.source}] {s.note}")
        needle = obs_map.get(key, m)
        try:
            status, ds = cdas.get_datasets(observatoryGroup=needle)
        except Exception as e:
            LOG.warning("live lookup failed for %s: %s", needle, e)
            continue
        print(f"    --- live CDAWeb datasets mentioning MAG/RPW/FIELDS for '{needle}' ---")
        for d in (ds or []):
            did = _first(d, "Id", "id") or ""
            label = _first(d, "Label", "label") or ""
            if any(t in (did + label).upper() for t in ("MAG", "RPW", "FLD", "FIELDS", "SCM")):
                print(f"    {did:34} {label[:60]}")
    return 0


# --------------------------------------------------------------------------- #
#  selftest -- exercise crash recovery without touching the network
# --------------------------------------------------------------------------- #
def selftest(_args) -> int:
    import tempfile
    LOG.info("self-test: manifest crash-recovery logic (no network)")
    d = tempfile.mkdtemp(prefix="sf_dl_test_")
    db = os.path.join(d, "manifest.db")
    m = C.ManifestDB(db)

    rows = [(f"http://x/file{i}.cdf", "DS", os.path.join(d, f"file{i}.cdf"), 100)
            for i in range(5)]
    m.upsert_pending(rows)
    assert len(m.get_work()) == 5, "all should start pending"

    # simulate: 2 finished (create their files), 1 interrupted mid-download
    for i in (0, 1):
        open(rows[i][2], "wb").write(b"x" * 100)
        m.mark(rows[i][0], "done", got=100)
    m.mark(rows[2][0], "downloading", bump_attempt=True)   # crashed here
    m.close()

    LOG.info("  --- simulating process restart ---")
    m2 = C.ManifestDB(db)
    rec = m2.reset_missing()                                # none missing yet
    work = m2.get_work()
    urls = {w[0] for w in work}
    ok = (len(work) == 3
          and rows[0][0] not in urls and rows[1][0] not in urls
          and rows[2][0] in urls)                           # interrupted -> retried
    LOG.info("  after restart: %d files queued (expect 3: 2 pending + 1 interrupted)",
             len(work))

    # simulate data loss: delete a 'done' file -> should be re-queued
    os.remove(rows[0][2])
    fixed = m2.reset_missing()
    requeued = rows[0][0] in {w[0] for w in m2.get_work()}
    LOG.info("  deleted a completed file -> reset_missing re-queued it: %s", requeued)
    m2.close()

    import shutil
    shutil.rmtree(d, ignore_errors=True)
    passed = ok and fixed == 1 and requeued
    LOG.info("PASS" if passed else "FAIL")
    return 0 if passed else 1


# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Download & sort EM CDF data for spacefoam analysis.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    f = sub.add_parser("fetch", help="discover + download data")
    f.add_argument("--missions", nargs="+",
                   help="any of: " + ", ".join(C.MISSION_REGISTRY))
    f.add_argument("--datasets", nargs="+",
                   help="restrict to specific CDAWeb dataset IDs")
    f.add_argument("--start", help="UTC start, e.g. 2022-01-01")
    f.add_argument("--end", help="UTC end, e.g. 2022-02-01")
    f.add_argument("--url-list", help="text file of CDF URLs (skips discovery)")
    f.add_argument("-o", "--output-dir", required=True, help="root output dir")
    f.add_argument("--workers", type=int, default=6, help="parallel downloads")
    f.add_argument("--chunk-bytes", type=int, default=1 << 16,
                   help="streaming chunk size (memory per worker)")
    f.add_argument("--timeout", type=float, default=60.0, help="per-request timeout (s)")
    f.add_argument("--dry-run", action="store_true", help="list files, download nothing")
    f.set_defaults(func=fetch)

    l = sub.add_parser("list-datasets", help="show current CDAWeb dataset IDs")
    l.add_argument("--missions", nargs="+", required=True)
    l.set_defaults(func=list_datasets)

    s = sub.add_parser("selftest", help="verify crash recovery (no network)")
    s.set_defaults(func=selftest)
    return p


def _install_signal_handlers():
    def handler(signum, frame):
        STOP.set()
    for sig in (signal.SIGINT, signal.SIGTERM):
        try:
            signal.signal(sig, handler)
        except Exception:
            pass


def main() -> int:
    _install_signal_handlers()
    args = build_parser().parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
