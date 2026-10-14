"""Builds ligint from the sources, for the test scripts when no prebuilt binary is given."""
import glob
import os
import subprocess

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))


def sources():
    return sorted(glob.glob(os.path.join(ROOT, "src", "core", "*.cpp"))) + [os.path.join(ROOT, "src", "cli", "main.cpp")]


def build_ligint(out, sanitize=False):
    """Builds the CLI into `out` with the compiler in $CXX (default g++). Returns any warnings."""
    cxx = os.environ.get("CXX", "g++")
    flags = (["-O1", "-g", "-fsanitize=address,undefined,float-cast-overflow", "-fno-sanitize-recover=all",
              "-fno-omit-frame-pointer"] if sanitize else ["-O2"])
    cmd = [cxx, "-std=c++17", "-Wall", "-Wextra", "-I", os.path.join(ROOT, "include"),
           "-I", os.path.join(ROOT, "src", "core")] + flags + sources() + ["-o", out]
    p = subprocess.run(cmd, capture_output=True)
    if p.returncode != 0:
        raise RuntimeError("build failed:\n" + p.stderr.decode())
    return p.stderr.decode().strip()
