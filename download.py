"""Download Voyager CDF products from CDAWeb in resumable time-chunks.

Each chunk is saved as netCDF under raw/<product>/<start>_<stop>.nc and skipped
on re-run, so an interrupted multi-decade pull picks up where it left off.
Live archive is NOT reachable from a sandbox; run this on a networked machine.
"""
from __future__ import annotations
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import List

from .config import Config, Product
from .catalog import resolve_dataset_ids, variables_for


def _chunks(t0: str, t1: str, days: int):
    fmt = "%Y-%m-%dT%H:%M:%SZ"
    a = datetime.strptime(t0, fmt).replace(tzinfo=timezone.utc)
    b = datetime.strptime(t1, fmt).replace(tzinfo=timezone.utc)
    step = timedelta(days=days)
    while a < b:
        c = min(a + step, b)
        yield a.strftime(fmt), c.strftime(fmt)
        a = c


def download(cfg: Config, only: List[str] | None = None, refresh_catalog: bool = False):
    from cdasws import CdasWs
    from cdasws.datarepresentation import DataRepresentation as DR
    cdas = CdasWs()
    cfg.paths.ensure()
    ids = resolve_dataset_ids(cfg, cdas, refresh=refresh_catalog)

    for p in cfg.products:
        if only and p.name not in only:
            continue
        ds_id = ids.get(p.name)
        if not ds_id:
            print(f"[{p.name}] no dataset id resolved; skipping.")
            continue
        outdir = cfg.paths.raw / p.name
        outdir.mkdir(parents=True, exist_ok=True)
        variables = variables_for(cdas, ds_id)
        print(f"[{p.name}] {ds_id}  ({len(variables)} vars)")

        for c0, c1 in _chunks(cfg.t_start, cfg.t_stop, cfg.download_chunk_days):
            fn = outdir / f"{c0[:10]}_{c1[:10]}.nc"
            if fn.exists():
                continue
            try:
                status, data = cdas.get_data(
                    ds_id, variables, c0, c1,
                    dataRepresentation=DR.XARRAY,
                )
            except Exception as e:
                print(f"   {c0[:10]} .. {c1[:10]}  ERROR {e}")
                continue
            if data is None or len(getattr(data, "dims", {})) == 0:
                continue  # genuine gap
            try:
                data.to_netcdf(fn)
                print(f"   {c0[:10]} .. {c1[:10]}  saved")
            except Exception as e:
                print(f"   {c0[:10]} .. {c1[:10]}  write-fail {e}")
