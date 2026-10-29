"""Download CRyPTIC source files over FTP with resume support and SHA-256 logging.

Files land in the raw directory exactly as served and are never rewritten once
complete. Partial transfers are kept as ``<name>.part`` and resumed with the FTP
REST command. Each completed file is recorded in ``download_log.json`` and
``SHA256SUMS`` alongside the data.
"""

from __future__ import annotations

import argparse
import datetime as dt
import ftplib
import hashlib
import json
import logging
from pathlib import Path
from urllib.parse import urljoin, urlparse

from tb_resistance.config import load_config, resolve

log = logging.getLogger(__name__)

LOG_NAME = "download_log.json"
BLOCK = 1 << 20


def sha256sum(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(BLOCK), b""):
            digest.update(block)
    return digest.hexdigest()


def _connect(host: str) -> ftplib.FTP:
    ftp = ftplib.FTP(host, timeout=120)
    ftp.login()
    ftp.voidcmd("TYPE I")
    return ftp


def fetch(url: str, dest: Path, retries: int = 5) -> Path:
    """Fetch ``url`` to ``dest``, resuming from ``dest.part`` if present."""
    if dest.exists():
        log.info("exists, skipping: %s", dest.name)
        return dest
    parsed = urlparse(url)
    part = dest.with_name(dest.name + ".part")
    dest.parent.mkdir(parents=True, exist_ok=True)

    for attempt in range(1, retries + 1):
        try:
            ftp = _connect(parsed.hostname)
            remote_size = ftp.size(parsed.path)
            offset = part.stat().st_size if part.exists() else 0
            if offset > remote_size:
                part.unlink()
                offset = 0
            if offset < remote_size:
                log.info("%s: %s/%s bytes, attempt %d", dest.name, offset, remote_size, attempt)
                with open(part, "ab") as fh:
                    ftp.retrbinary(f"RETR {parsed.path}", fh.write, blocksize=BLOCK,
                                   rest=offset or None)
            ftp.close()
            if part.stat().st_size != remote_size:
                raise OSError(f"size mismatch for {dest.name}")
            part.rename(dest)
            return dest
        except (ftplib.all_errors) as exc:  # network errors: retry from where we stopped
            log.warning("%s failed (%s); retrying", dest.name, exc)
    raise RuntimeError(f"could not download {url} after {retries} attempts")


def _load_log(raw_dir: Path) -> dict:
    path = raw_dir / LOG_NAME
    return json.loads(path.read_text()) if path.exists() else {}


def _write_log(raw_dir: Path, entries: dict) -> None:
    (raw_dir / LOG_NAME).write_text(json.dumps(entries, indent=2) + "\n")
    lines = [f"{e['sha256']}  {name}" for name, e in sorted(entries.items())]
    (raw_dir / "SHA256SUMS").write_text("\n".join(lines) + "\n")


def download_all(cfg: dict) -> dict:
    """Download every file listed in ``source.files`` and return the log entries."""
    base = cfg["source"]["base_url"]
    raw_dir = resolve(cfg["paths"]["raw_dir"])
    entries = _load_log(raw_dir)
    for key, rel in cfg["source"]["files"].items():
        url = urljoin(base, rel)
        dest = raw_dir / Path(rel).name
        fetch(url, dest)
        if dest.name not in entries:
            entries[dest.name] = {
                "key": key,
                "url": url,
                "bytes": dest.stat().st_size,
                "sha256": sha256sum(dest),
                "downloaded_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
            }
            _write_log(raw_dir, entries)
            log.info("%s sha256=%s", dest.name, entries[dest.name]["sha256"])
    return entries


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/data.yaml")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    download_all(load_config(resolve(args.config)))


if __name__ == "__main__":
    main()
