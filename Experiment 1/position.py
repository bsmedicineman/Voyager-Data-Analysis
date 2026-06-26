"""Where in space each anomaly happened, via SPICE.

Heliocentric distance (AU) for Voyager 1 (NAIF -31) / Voyager 2 (-32) relative
to the Sun. Needs three kernels, fetched once into kernels/:
  * naif0012.tls            leap seconds
  * de440.bsp               planetary + Sun ephemeris
  * voyager_1.ST+....bsp    spacecraft trajectory (likewise voyager_2)

Kernel URLs live at https://naif.jpl.nasa.gov/pub/naif/  (not reachable from a
sandbox; fetch on a networked machine). Exact Voyager SPK filenames change, so
`fetch_kernels` lists the directory and grabs the newest matching files.
"""
from __future__ import annotations
from pathlib import Path
from typing import Dict, List
import numpy as np

NAIF_ID = {"voyager1": -31, "voyager2": -32}
NAIF_BASE = "https://naif.jpl.nasa.gov/pub/naif/"
GENERIC = NAIF_BASE + "generic_kernels/"
VGR = NAIF_BASE + "VOYAGER/kernels/spk/"


def fetch_kernels(kernel_dir: Path) -> None:
    """Download leap-second, planetary, and both Voyager SPK kernels."""
    import urllib.request, re
    kernel_dir.mkdir(parents=True, exist_ok=True)

    def get(url, name):
        dst = kernel_dir / name
        if dst.exists():
            return
        print(f"  fetching {name}")
        urllib.request.urlretrieve(url, dst)

    get(GENERIC + "lsk/naif0012.tls", "naif0012.tls")
    get(GENERIC + "spk/planets/de440.bsp", "de440.bsp")
    # Voyager SPKs: list dir, take newest *.bsp per spacecraft
    html = urllib.request.urlopen(VGR).read().decode("utf-8", "ignore")
    files = re.findall(r'href="([^"]+\.bsp)"', html)
    for sc, tag in (("voyager_1", "voyager1"), ("voyager_2", "voyager2")):
        cands = sorted(f for f in files if sc in f.lower())
        if cands:
            get(VGR + cands[-1], cands[-1])


def furnish(kernel_dir: Path) -> None:
    import spiceypy as spice
    for k in sorted(kernel_dir.glob("*")):
        if k.suffix.lower() in (".tls", ".bsp", ".tpc", ".tf", ".tsc"):
            spice.furnsh(str(k))


def heliocentric_au(epoch_seconds: np.ndarray, spacecraft: str,
                    kernel_dir: Path) -> np.ndarray:
    """AU distance from Sun at each epoch-second time. NaN where outside coverage."""
    import spiceypy as spice
    furnish(kernel_dir)
    targ = str(NAIF_ID[spacecraft])
    out = np.full(len(epoch_seconds), np.nan)
    # epoch seconds (since 1970) -> ET
    for i, es in enumerate(np.atleast_1d(epoch_seconds)):
        try:
            utc = np.datetime64(int(es), "s").astype("datetime64[ms]").astype(str)
            et = spice.str2et(utc.replace("T", " "))
            pos, _ = spice.spkpos(targ, et, "ECLIPJ2000", "NONE", "SUN")
            out[i] = spice.convrt(float(np.linalg.norm(pos)), "km", "AU")
        except Exception:
            pass
    return out


def attach_distance(rows: List[Dict], kernel_dir: Path) -> List[Dict]:
    """Add 'distance_au' to each anomaly row (uses its midpoint time)."""
    by_sc: Dict[str, List[int]] = {}
    for i, r in enumerate(rows):
        by_sc.setdefault(r["spacecraft"], []).append(i)
    for sc, idx in by_sc.items():
        mid = np.array([(rows[i]["t_start"] + rows[i]["t_stop"]) / 2 for i in idx])
        d = heliocentric_au(mid, sc, kernel_dir)
        for k, i in enumerate(idx):
            rows[i]["distance_au"] = float(d[k]) if np.isfinite(d[k]) else None
    return rows
