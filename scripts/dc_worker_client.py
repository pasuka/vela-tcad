"""Sequential, bounded JSON-lines transport for vela_example_runner --dc-worker."""
import json
import queue
import subprocess
import threading


class DCWorker:
    def __init__(self, runner, cwd, env, cpu_reader, *, reuse_linear_analysis=True):
        self.log = (cwd / 'worker_stderr.log').open('w', encoding='utf-8')
        self.proc = subprocess.Popen(
            [str(runner), '--dc-worker' if reuse_linear_analysis else '--dc-worker-no-linear-reuse'], cwd=cwd, env=env,
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=self.log,
            text=True, encoding='utf-8', bufsize=1)
        self.cpu_reader = cpu_reader
        self.count = 0
        self.responses = queue.Queue()
        def collect():
            try:
                for line in self.proc.stdout:
                    self.responses.put(line)
            finally:
                self.responses.put(None)
        self.reader = threading.Thread(target=collect, daemon=True)
        self.reader.start()

    def run(self, config, timeout=600):
        before = self.cpu_reader(self.proc) if self.count else dict(kernel=0., user=0., total=0.)
        self.count += 1
        self.proc.stdin.write(json.dumps({'id': self.count, 'config': str(config)}) + '\n')
        self.proc.stdin.flush()
        try:
            line = self.responses.get(timeout=timeout)
        except queue.Empty as error:
            self.proc.kill()
            self.proc.wait()
            raise TimeoutError('DC worker request timed out') from error
        if line is None:
            raise RuntimeError(f'DC worker exited before response {self.count}')
        response = json.loads(line)
        if response.get('id') != self.count:
            raise RuntimeError('DC worker response ID mismatch')
        after = self.cpu_reader(self.proc)
        cpu = {k: after[k] - before[k] for k in before}
        return response, cpu

    def close(self):
        if self.proc.poll() is None:
            try:
                self.proc.stdin.write('{"shutdown":true}\n')
                self.proc.stdin.flush()
                self.proc.wait(timeout=10)
            except (OSError, subprocess.TimeoutExpired):
                self.proc.kill()
                self.proc.wait()
        self.proc.stdin.close()
        self.reader.join(timeout=10)
        self.proc.stdout.close()
        self.log.close()
