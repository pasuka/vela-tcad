"""Working-tree snapshot semantics used by Invoke-VelaVm.ps1."""
import pathlib
import subprocess
import sys
import tarfile
import tempfile
import unittest

SCRIPT = pathlib.Path(__file__).resolve().parents[2] / "scripts/remote-vm/snapshot.py"


class SnapshotTests(unittest.TestCase):
    def test_local_edits_new_files_deletions_and_ignored_outputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary) / "repo"
            root.mkdir()
            subprocess.run(["git", "init", "-q", str(root)], check=True)
            (root / ".gitignore").write_text("build/\n", encoding="utf-8")
            (root / "tracked.txt").write_text("original", encoding="utf-8")
            (root / "deleted.txt").write_text("deleted", encoding="utf-8")
            subprocess.run(["git", "add", "."], cwd=root, check=True)
            (root / "tracked.txt").write_text("working edit", encoding="utf-8")
            (root / "deleted.txt").unlink()
            (root / "new file.txt").write_text("new content", encoding="utf-8")
            (root / "build").mkdir()
            (root / "build/output.txt").write_text("ignored", encoding="utf-8")
            archive = pathlib.Path(temporary) / "snapshot.tar.gz"
            subprocess.run([sys.executable, str(SCRIPT), str(root), str(archive)], check=True)
            with tarfile.open(archive) as package:
                self.assertEqual(set(package.getnames()), {".gitignore", "tracked.txt", "new file.txt"})
                self.assertEqual(package.extractfile("tracked.txt").read(), b"working edit")
                self.assertEqual(package.extractfile("new file.txt").read(), b"new content")


if __name__ == "__main__":
    unittest.main()
