"""Archive the current working tree's Git-visible files, including local edits."""
import pathlib
import subprocess
import sys
import tarfile


def main():
    root = pathlib.Path(sys.argv[1]).resolve()
    destination = pathlib.Path(sys.argv[2]).resolve()
    names = subprocess.check_output(
        ["git", "-c", "core.fsmonitor=false", "ls-files", "-z", "--cached",
         "--others", "--exclude-standard"], cwd=root
    ).decode("utf-8").split("\0")
    count = 0
    with tarfile.open(destination, "w:gz") as archive:
        for name in sorted(set(filter(None, names))):
            path = root / name
            if not path.exists() and not path.is_symlink():
                continue  # Deleted tracked files must stay deleted in the snapshot.
            if path.is_symlink() or path.is_dir():
                raise RuntimeError(f"Snapshot requires regular files; inspect {name}")
            archive.add(path, arcname=name, recursive=False)
            count += 1
    print(f"Snapshot: {count} files, {destination.stat().st_size} bytes")


if __name__ == "__main__":
    main()
