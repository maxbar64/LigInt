"""Compilazione di ligint dai sorgenti, per gli script di test quando non si passa un binario gia' pronto."""
import glob
import os
import subprocess

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))


def sources():
    return sorted(glob.glob(os.path.join(ROOT, "src", "core", "*.cpp"))) + [os.path.join(ROOT, "src", "cli", "main.cpp")]


def build_ligint(out, sanitize=False):
    """Compila la CLI in `out` con il compilatore in $CXX (default g++). Ritorna gli eventuali avvisi."""
    cxx = os.environ.get("CXX", "g++")
    flags = (["-O1", "-g", "-fsanitize=address,undefined,float-cast-overflow", "-fno-sanitize-recover=all",
              "-fno-omit-frame-pointer"] if sanitize else ["-O2"])
    cmd = [cxx, "-std=c++17", "-Wall", "-Wextra", "-I", os.path.join(ROOT, "include"),
           "-I", os.path.join(ROOT, "src", "core")] + flags + sources() + ["-o", out]
    p = subprocess.run(cmd, capture_output=True)
    if p.returncode != 0:
        raise RuntimeError("compilazione fallita:\n" + p.stderr.decode())
    return p.stderr.decode().strip()
