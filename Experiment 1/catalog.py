"""Resolve CDAWeb dataset IDs from instrument keywords.

Exact CDAWeb ids change over time, so rather than hard-code possibly-stale
strings we query the live catalog and match on keywords. Results are cached
to disk so discovery only hits the network once.
"""
from __future__ import annotations
import json
from typing import Dict, List, Optional

from .config import Config, Product


def _all_voyager_datasets(cdas) -> List[Dict]:
    out = []
    for grp in ("Voyager",):
        try:
            out += cdas.get_datasets(observatoryGroup=grp)
        except Exception:
            pass
    if not out:  # fallback: pull everything and filter by id prefix
        try:
            out = [d for d in cdas.get_datasets()
                   if "VOYAGER" in d.get("Id", "").upper()]
        except Exception:
            out = []
    return out


def _score(dataset: Dict, keywords: List[str]) -> int:
    hay = (dataset.get("Id", "") + " " + dataset.get("Label", "")).lower()
    return sum(1 for k in keywords if k.lower() in hay)


def resolve_dataset_ids(cfg: Config, cdas, refresh: bool = False) -> Dict[str, str]:
    """Return {product.name: cdaweb_dataset_id}. Caches to tables/datasets.json."""
    cfg.paths.ensure()
    cache = cfg.paths.tables / "datasets.json"
    resolved: Dict[str, str] = {}
    if cache.exists() and not refresh:
        resolved = json.loads(cache.read_text())

    catalog = None
    for p in cfg.products:
        if p.name in resolved and p.dataset_id is None:
            continue
        if p.dataset_id:
            resolved[p.name] = p.dataset_id
            continue
        if catalog is None:
            catalog = _all_voyager_datasets(cdas)
        ranked = sorted(catalog, key=lambda d: _score(d, p.catalog_keywords),
                        reverse=True)
        best = ranked[0] if ranked and _score(ranked[0], p.catalog_keywords) >= 2 else None
        if best:
            resolved[p.name] = best["Id"]
        else:
            print(f"  [warn] could not auto-resolve dataset for {p.name}; "
                  f"set Product.dataset_id manually. keywords={p.catalog_keywords}")

    cache.write_text(json.dumps(resolved, indent=2))
    return resolved


def variables_for(cdas, dataset_id: str) -> List[str]:
    """All variable names in a dataset (we let analysis pick the relevant ones)."""
    try:
        return [v["Name"] for v in cdas.get_variables(dataset_id)]
    except Exception:
        return ["ALL-VARIABLES"]
