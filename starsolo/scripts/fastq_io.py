"""Minimal paired-FASTQ I/O over (pi)gzip pipes -- no compiled dependencies.

Replaces dnaio, whose compression libraries could not be installed on the
cluster. Decompression and compression run in external pigz (or gzip)
processes, so the Python side only splits lines.

Records are bytes, each line WITH its trailing newline.
"""
import os
import shutil
import subprocess

BUF = 1 << 20


def _gzip_tool():
    tool = shutil.which("pigz") or shutil.which("gzip")
    if not tool:
        raise RuntimeError("neither pigz nor gzip is on PATH")
    return tool


class _Reader:
    def __init__(self, path):
        self.path = path
        if path.endswith(".gz"):
            self.proc = subprocess.Popen([_gzip_tool(), "-dc", path],
                                         stdout=subprocess.PIPE, bufsize=BUF)
            self.fh = self.proc.stdout
        else:
            self.proc = None
            self.fh = open(path, "rb", buffering=BUF)

    def close(self, complete):
        if self.proc and not complete:
            self.proc.kill()            # stopped early on purpose
        self.fh.close()
        if self.proc:
            rc = self.proc.wait()
            if complete and rc != 0:
                raise IOError(f"decompressing {self.path} failed (exit {rc})")


def _name(header):
    name = header.split(None, 1)[0]
    return name[:-2] if name.endswith((b"/1", b"/2")) else name


def read_pairs(path1, path2, check_every=100_000):
    """Yield (h1, s1, q1, h2, s2, q2) for each read pair.

    Fails if the files have different read counts, a record is truncated, or
    read names disagree (checked on the first pair and every `check_every`).
    Wrap in contextlib.closing() when not consuming everything.
    """
    r1, r2 = _Reader(path1), _Reader(path2)
    complete = False
    try:
        l1, l2 = r1.fh.readline, r2.fh.readline
        n = 0
        while True:
            h1 = l1()
            if not h1:
                if l2():
                    raise ValueError(f"{path2} has more reads than {path1}")
                complete = True
                return
            s1 = l1(); l1(); q1 = l1()
            h2 = l2(); s2 = l2(); l2(); q2 = l2()
            if not q1 or not q2:
                raise ValueError(f"truncated record or unequal read counts at pair {n + 1}")
            if n % check_every == 0 and _name(h1) != _name(h2):
                raise ValueError(f"read names differ at pair {n + 1}: {h1!r} vs {h2!r}")
            n += 1
            yield h1, s1, q1, h2, s2, q2
    finally:
        r1.close(complete)
        r2.close(complete)


class Writer:
    """Write bytes to a file, through `pigz -p threads` / gzip when it ends in .gz."""

    def __init__(self, path, threads=1, level=1):
        self.path = path
        self.proc = None
        if path.endswith(".gz"):
            tool = _gzip_tool()
            cmd = [tool, f"-{level}", "-c"]
            if os.path.basename(tool) == "pigz":
                cmd += ["-p", str(max(1, threads))]
            self._out = open(path, "wb")
            self.proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=self._out,
                                         bufsize=BUF)
            self.write = self.proc.stdin.write
        else:
            self._out = open(path, "wb", buffering=BUF)
            self.write = self._out.write

    def close(self):
        if self.proc:
            self.proc.stdin.close()
            rc = self.proc.wait()
            self._out.close()
            if rc != 0:
                raise IOError(f"compressing {self.path} failed (exit {rc})")
        else:
            self._out.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
