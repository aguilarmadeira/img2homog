"""Minimal tests: run with  python -m pytest tests/  (or just python tests/test_img2homog.py)."""
import subprocess, sys, os, numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def run(args):
    return subprocess.run([sys.executable, os.path.join(HERE, "img2homog.py")] + args,
                          capture_output=True, text=True, cwd=HERE)

def test_examples_pass():
    ref = {"reentrant_v1": 0.398, "missing_rib_v1": 0.345,
           "chiral_v1": 0.289, "rotating_v1": 0.397}
    for name, frac in ref.items():
        r = run([f"examples/{name}.png", "-N", "100",
                 "--reference-fraction", str(frac), "-o", f"/tmp/{name}.txt"])
        assert r.returncode == 0, r.stdout + r.stderr
        A = np.loadtxt(f"/tmp/{name}.txt")
        assert A.shape == (100, 100)
        assert abs(A.mean() - frac) < 0.02, (name, A.mean())
        assert "[FAIL]" not in r.stdout

def test_disconnected_fails():
    from PIL import Image
    A = np.ones((50, 50), np.uint8) * 255
    A[10:20, 10:20] = 0   # isolated solid island, no periodic connection
    Image.fromarray(A).save("/tmp/island.png")
    r = run(["/tmp/island.png", "-N", "50", "-o", "/tmp/island.txt"])
    assert r.returncode == 1
    assert "[FAIL]" in r.stdout

if __name__ == "__main__":
    test_examples_pass(); test_disconnected_fails(); print("all tests passed")
