import tarfile
import logging
from urllib.parse import urlsplit
from urllib.request import urlopen, Request
from pathlib import Path


def download_and_unpack_or_fetch(
    url: str,
    dest_dir: str,
    logger: logging.Logger
) -> str:
    """
    Download from `url` into `dest_dir` (creating it if needed).

    - Archives (.tar.gz/.tgz):
        • If already downloaded & unpacked, returns existing unpacked path.
        • Otherwise downloads & extracts, then returns unpacked folder path.
    - Plain files:
        • If already downloaded, returns existing file path.
        • Otherwise downloads it, then returns the file path.

    Logs each step and shows download progress via the given logger.

    :param url: HTTP(S) URL of the file or archive to download.
    :param dest_dir: Directory to save the downloaded file or unpacked archive.
    :param logger: Logger for status/progress messages.
    :return: Absolute path to the downloaded file or unpacked archive.
    """
    dest = Path(dest_dir).expanduser().absolute()
    logger.info(f"Using destination directory: {dest}")
    dest.mkdir(parents=True, exist_ok=True)

    filename = Path(urlsplit(url).path).name
    local_path = dest / filename
    logger.info(f"Local path will be: {local_path}")

    def _download():
        req = Request(url, headers={'User-Agent': 'Python urllib'})
        with urlopen(req) as resp, open(local_path, 'wb') as out:
            total_size = int(resp.getheader('Content-Length') or 0)
            if total_size:
                logger.info(f"Starting download of {filename} ({total_size} bytes)")
                next_pct = 10
            else:
                logger.info(f"Starting download of {filename} (size unknown)")
                threshold = 1 * 1024 * 1024  # 1 MiB
                next_bytes = threshold

            downloaded = 0
            chunk_size = 8192

            while True:
                chunk = resp.read(chunk_size)
                if not chunk:
                    break
                out.write(chunk)
                downloaded += len(chunk)

                if total_size:
                    pct = downloaded * 100 / total_size
                    if pct >= next_pct:
                        logger.info(f"Downloaded {downloaded}/{total_size} bytes ({pct:.0f}%)")
                        next_pct += 10
                else:
                    if downloaded >= next_bytes:
                        logger.info(f"Downloaded {downloaded} bytes")
                        next_bytes += threshold

            logger.info("Download complete")

    is_tar = filename.lower().endswith(('.tar.gz', '.tgz'))

    if is_tar:
        logger.info(f"Detected archive format for {filename}")

        def _get_roots(tar_path):
            with tarfile.open(tar_path, 'r:gz') as tar:
                return {Path(m.name).parts[0] for m in tar.getmembers() if m.name.strip()}

        # If archive exists, see if its contents are already unpacked
        if local_path.exists():
            logger.info(f"Archive already exists at {local_path}, checking extraction")
            try:
                roots = _get_roots(local_path)
                if len(roots) == 1:
                    candidate = dest / roots.pop()
                    if candidate.exists():
                        logger.info(f"Found existing extracted folder {candidate}, skipping.")
                        return str(candidate)
                else:
                    others = [p for p in dest.iterdir() if p != local_path]
                    if others:
                        logger.info(f"Found already-extracted contents in {dest}, skipping.")
                        return str(dest)
            except tarfile.TarError:
                logger.warning(f"Archive at {local_path} seems corrupt. Re-downloading.")

        # Download & extract
        _download()
        logger.info(f"Extracting {local_path} to {dest}")
        with tarfile.open(local_path, 'r:gz') as tar:
            roots = {Path(m.name).parts[0] for m in tar.getmembers() if m.name.strip()}
            tar.extractall(path=dest)

        if len(roots) == 1:
            unpacked = dest / roots.pop()
            logger.info(f"Extraction complete: single folder {unpacked}")
            return str(unpacked)
        else:
            logger.info(f"Extraction complete: multiple items, returning {dest}")
            return str(dest)

    else:
        # Plain file
        logger.info(f"Detected single file for {filename}")
        if not local_path.exists():
            logger.info(f"{local_path} not found locally; downloading.")
            _download()
        else:
            logger.info(f"File already exists at {local_path}, skipping download.")
        return str(local_path)
