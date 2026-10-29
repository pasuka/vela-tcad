"""Resource paths in a source checkout or a self-contained zipapp."""
from __future__ import annotations
import atexit
from importlib.resources import files
from pathlib import Path
import tempfile

_temporary = None


def resource_root() -> Path:
    root = Path(__file__).resolve().parents[2]
    if (root / "schemas").is_dir():
        return root
    global _temporary
    if _temporary is None:
        _temporary = tempfile.TemporaryDirectory(prefix="vela_sdevice_resources_")
        atexit.register(_temporary.cleanup)
        def copy_tree(source, target):
            target.mkdir(parents=True, exist_ok=True)
            for item in source.iterdir():
                path = target / item.name
                if item.is_dir():
                    copy_tree(item, path)
                else:
                    path.write_bytes(item.read_bytes())
        copy_tree(files("vela_sdevice").joinpath("_assets"), Path(_temporary.name))
    return Path(_temporary.name)
