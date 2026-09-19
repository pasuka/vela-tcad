"""Read-only Windows aggregate CPU sampling, without changing other processes."""
import ctypes
from ctypes import wintypes
import time


def cpu_snapshot():
    idle,kernel,user=(wintypes.FILETIME() for _ in range(3))
    api=ctypes.WinDLL('kernel32',use_last_error=True).GetSystemTimes
    api.argtypes=[ctypes.POINTER(wintypes.FILETIME)]*3
    api.restype=wintypes.BOOL
    if not api(ctypes.byref(idle),ctypes.byref(kernel),ctypes.byref(user)):
        raise ctypes.WinError(ctypes.get_last_error())
    ticks=lambda t:(t.dwHighDateTime<<32)|t.dwLowDateTime
    return dict(monotonic=time.perf_counter(),idle=ticks(idle),total=ticks(kernel)+ticks(user))


def cpu_interval(a,b):
    total=b['total']-a['total']; idle=b['idle']-a['idle']
    if total<=0 or idle<0 or idle>total:raise ValueError('Invalid CPU sampling interval')
    return dict(interval_seconds=b['monotonic']-a['monotonic'],
                busy_percent=100.*(total-idle)/total,
                busy_cpu_seconds=(total-idle)/1e7)


def wait_for_idle(limit,seconds=15,timeout=600):
    """Require three consecutive five-second windows; preserve all observations."""
    start=time.perf_counter();previous=cpu_snapshot();quiet=0.;samples=[]
    while quiet<seconds:
        if time.perf_counter()-start>timeout:
            error=TimeoutError('CPU idle gate timed out');error.samples=samples;raise error
        time.sleep(5)
        current=cpu_snapshot();sample=cpu_interval(previous,current);samples.append(sample)
        quiet=quiet+sample['interval_seconds'] if sample['busy_percent']<=limit else 0.
        previous=current
        if len(samples)%6==0:
            print(f'CPU gate: last={sample["busy_percent"]:.1f}%, required<={limit:.1f}%, quiet={quiet:.1f}s',flush=True)
    return samples
