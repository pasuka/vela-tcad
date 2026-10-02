"""Lock ownership, atomic records and numbered frontend failures."""
from __future__ import annotations
import ctypes
import json
import os
from pathlib import Path
import socket
import time
import uuid


class RunError(Exception):
    def __init__(self, message_id, message, **details):
        self.message_id, self.message, self.details = message_id, message, details
        self.exit_code = int(message_id[1])
        super().__init__(f"{message_id}: {message}")


def write_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + ".tmp")
    try:
        temporary.write_text(json.dumps(value, indent=2, ensure_ascii=True, allow_nan=False) + "\n", encoding="utf-8")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def pid_alive(pid):
    if not isinstance(pid, int) or pid <= 0:
        return True
    if os.name == "nt":
        from ctypes import wintypes
        kernel = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
        kernel.OpenProcess.restype = wintypes.HANDLE
        kernel.GetExitCodeProcess.argtypes = [wintypes.HANDLE, ctypes.POINTER(wintypes.DWORD)]
        kernel.CloseHandle.argtypes = [wintypes.HANDLE]
        handle = kernel.OpenProcess(0x1000, False, pid)
        if not handle:
            return ctypes.get_last_error() != 87  # invalid PID; access denied stays live
        try:
            code = wintypes.DWORD()
            return not kernel.GetExitCodeProcess(handle, ctypes.byref(code)) or code.value == 259
        finally:
            kernel.CloseHandle(handle)
    try:
        os.kill(pid, 0)
        return True
    except ProcessLookupError:
        return False
    except PermissionError:
        return True


class RunLock:
    def __init__(self, path):
        self.path = Path(path)
        self.owner = {"pid": os.getpid(), "host": socket.gethostname(), "time": time.time(), "token": uuid.uuid4().hex}
        self.owned = False
        self.stale = False
        self.previous = None

    def acquire(self):
        try:
            self._create()
        except FileExistsError:
            try:
                previous = json.loads(self.path.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                raise RunError("E701", f"Another vela device run holds '{self.path}' (unreadable ownership)")
            if previous.get("host") != self.owner["host"] or pid_alive(previous.get("pid")):
                raise RunError("E701", f"Another vela device run holds '{self.path}'", owner=previous)
            # One reclaimer owns this guard, preventing two stale-lock readers
            # from unlinking each other's new live lock.
            guard = self.path.with_name(".lock.reclaim")
            try:
                fd = os.open(guard, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
            except FileExistsError:
                raise RunError("E701", f"Another run is reclaiming '{self.path}'")
            try:
                os.close(fd)
                if json.loads(self.path.read_text(encoding="utf-8")) != previous:
                    raise RunError("E701", f"Lock ownership changed for '{self.path}'")
                self.path.unlink()
                try:
                    self._create()
                except FileExistsError:
                    raise RunError("E701", f"Another vela device run holds '{self.path}'")
                self.stale = True
                self.previous = previous
            finally:
                guard.unlink(missing_ok=True)
        return self

    def _create(self):
        fd = os.open(self.path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            json.dump(self.owner, stream)
            stream.flush()
            os.fsync(stream.fileno())
        self.owned = True

    def release(self):
        if self.owned:
            try:
                if json.loads(self.path.read_text(encoding="utf-8")) == self.owner:
                    self.path.unlink()
            except FileNotFoundError:
                pass
            self.owned = False
