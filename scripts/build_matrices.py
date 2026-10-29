"""Build per-drug feature matrices and labels from the CRyPTIC release.

Usage:
    python scripts/build_matrices.py [--config configs/data.yaml] [--skip-download]
                                     [--rebuild-interim]
"""

from __future__ import annotations

import argparse
import json
import logging

from tb_resistance.build import build_all
from tb_resistance.config import load_config, resolve
from tb_resistance.download import download_all


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default="configs/data.yaml")
    parser.add_argument("--skip-download", action="store_true",
                        help="use files already in the raw dir")
    parser.add_argument("--rebuild-interim", action="store_true",
                        help="re-filter the raw MUTATIONS table even if a cache exists")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    cfg = load_config(resolve(args.config))
    if not args.skip_download:
        download_all(cfg)
    report = build_all(cfg, rebuild_interim=args.rebuild_interim)
    print(json.dumps(report, indent=2, default=str))


if __name__ == "__main__":
    main()
