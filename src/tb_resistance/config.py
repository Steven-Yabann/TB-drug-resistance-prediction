"""Load the YAML data configuration and resolve repo-relative paths."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = REPO_ROOT / "configs" / "data.yaml"


def load_config(path: str | Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    for key, value in cfg["paths"].items():
        if Path(value).is_absolute():
            raise ValueError(f"paths.{key} must be repo-relative, got {value!r}")
    return cfg


def resolve(relative: str | Path) -> Path:
    """Resolve a repo-relative path against the repository root."""
    return REPO_ROOT / relative
