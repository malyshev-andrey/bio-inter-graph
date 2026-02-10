import json
import threading
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

from biointergraph.shared import remote_file2local


class _QuietHandler(SimpleHTTPRequestHandler):
    def log_message(self, format, *args):  # pragma: no cover - keep test output clean
        pass


def _start_http_server(directory: Path):
    handler = partial(_QuietHandler, directory=str(directory))
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    base_url = f"http://127.0.0.1:{server.server_port}"
    return server, thread, base_url


def test_remote_file2local_downloads_and_caches(tmp_path: Path) -> None:
    remote_dir = tmp_path / "remote"
    remote_dir.mkdir()
    (remote_dir / "data.txt").write_text("hello", encoding="utf-8")

    server, _thread, base_url = _start_http_server(remote_dir)
    try:
        url = f"{base_url}/data.txt"
        cache_dir = tmp_path / "cache"

        new_url = remote_file2local(url, cache_dir=str(cache_dir))
    finally:
        server.shutdown()
        server.server_close()

    assert new_url.startswith("file://")
    local_path = Path(new_url.replace("file://", "", 1))
    assert local_path.read_text(encoding="utf-8") == "hello"
    assert local_path.is_file()

    metadata_path = local_path.with_suffix(local_path.suffix + ".meta.json")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    assert metadata["url"] == url
    assert metadata["local_path"] == str(local_path)
    assert metadata["download_time"] >= 0


def test_remote_file2local_uses_cache_without_remote(tmp_path: Path) -> None:
    remote_dir = tmp_path / "remote"
    remote_dir.mkdir()
    (remote_dir / "data.txt").write_text("cached", encoding="utf-8")

    server, _thread, base_url = _start_http_server(remote_dir)
    url = f"{base_url}/data.txt"
    cache_dir = tmp_path / "cache"

    first_url = remote_file2local(url, cache_dir=str(cache_dir))
    local_path = Path(first_url.replace("file://", "", 1))

    server.shutdown()
    server.server_close()
    (remote_dir / "data.txt").unlink()

    second_url = remote_file2local(url, cache_dir=str(cache_dir))

    assert second_url == first_url
    assert local_path.read_text(encoding="utf-8") == "cached"


def test_remote_file2local_preserves_protocol_chain(tmp_path: Path) -> None:
    remote_dir = tmp_path / "remote"
    remote_dir.mkdir()
    (remote_dir / "data.txt").write_text("chain", encoding="utf-8")

    server, _thread, base_url = _start_http_server(remote_dir)
    try:
        url = f"{base_url}/data.txt"
        cache_dir = tmp_path / "cache"

        chained_url = f"simplecache::{url}"
        new_url = remote_file2local(chained_url, cache_dir=str(cache_dir))
    finally:
        server.shutdown()
        server.server_close()

    assert new_url.startswith("simplecache::file://")
    local_path = Path(new_url.replace("simplecache::file://", "", 1))
    assert local_path.read_text(encoding="utf-8") == "chain"
