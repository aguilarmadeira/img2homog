#!/usr/bin/env python3
"""img2homog — black-and-white image to 0/1 periodic bitmap for
finite-element asymptotic homogenization (Guedes & Kikuchi, 1990).

Converts a raster image of a periodic unit cell into an N x N matrix of
0/1 values (1 = solid, 0 = void) and runs the validation checks that
decide whether the bitmap is actually usable for homogenization:

  * solid fraction (optionally compared with a reference value);
  * periodic 4-connectivity: a single solid component that spans the
    cell boundaries in both directions;
  * isolated pixels and diagonal-only contacts (spurious point hinges);
  * thinnest-member estimate, with a warning below 4 pixels;
  * JPEG input warning (compression noise near edges).

Orientation convention (IMPORTANT for chiral cells):
  row 1 of the output = TOP of the image (y max); column 1 = LEFT (x min).
  Flipping rows reverses chirality and silently changes the sign of the
  normal-shear coupling of chiral lattices. Keep the convention consistent
  between this tool and the homogenization code.

Usage:
  python img2homog.py cell.png -N 200 -o cell_200.txt
  python img2homog.py cell.png -N 100 --solid dark --threshold 0.5 \
      --reference-fraction 0.345 --preview cell_100.png

Author: J. F. A. Madeira (IDMEC, Instituto Superior Tecnico / ISEL).
License: MIT.
"""
import argparse, os, sys
import numpy as np
from PIL import Image
from scipy import ndimage


def load_binary(path, solid, threshold):
    img = Image.open(path).convert("L")
    A = np.asarray(img, dtype=float) / 255.0
    warn_jpeg = os.path.splitext(path)[1].lower() in (".jpg", ".jpeg")
    if solid == "dark":
        B = A < threshold
    else:
        B = A >= threshold
    grey = np.mean((A > 0.02) & (A < 0.98))
    return B.astype(np.uint8), warn_jpeg, grey


def resample(B, N, rule, min_frac):
    """Resample binary image to N x N. rule='area': pixel=1 iff solid area
    fraction inside the target pixel >= min_frac; rule='center': pixel=1 iff
    the source pixel containing the target-pixel centre is solid."""
    H, W = B.shape
    out = np.zeros((N, N), dtype=np.uint8)
    if rule == "center":
        for i in range(N):
            for j in range(N):
                y = int((i + 0.5) * H / N)
                x = int((j + 0.5) * W / N)
                out[i, j] = B[min(y, H-1), min(x, W-1)]
        return out
    # area rule via integral image
    S = np.zeros((H + 1, W + 1))
    S[1:, 1:] = np.cumsum(np.cumsum(B, 0), 1)
    ys = np.linspace(0, H, N + 1)
    xs = np.linspace(0, W, N + 1)
    for i in range(N):
        y0, y1 = ys[i], ys[i + 1]
        iy0, iy1 = int(np.floor(y0)), int(np.ceil(y1))
        for j in range(N):
            x0, x1 = xs[j], xs[j + 1]
            ix0, ix1 = int(np.floor(x0)), int(np.ceil(x1))
            # integer-box approximation of the solid area fraction
            area = S[iy1, ix1] - S[iy0, ix1] - S[iy1, ix0] + S[iy0, ix0]
            out[i, j] = 1 if area / max((iy1 - iy0) * (ix1 - ix0), 1) >= min_frac else 0
    return out


def checks(A):
    N = A.shape[0]
    rep = {}
    rep["solid_fraction"] = float(A.mean())
    T = np.tile(A, (3, 3))
    lab, n = ndimage.label(T)                      # 4-connectivity
    centre = lab[N:2 * N, N:2 * N]
    ids = sorted(set(centre[centre > 0].ravel()))
    rep["single_component"] = len(ids) == 1
    rep["n_components_in_cell"] = len(ids)
    rep["periodic_span_x"] = rep["periodic_span_y"] = False
    if ids:
        mask = (lab == ids[0])
        rep["periodic_span_y"] = bool(mask[:N, :].any() and mask[2 * N:, :].any())
        rep["periodic_span_x"] = bool(mask[:, :N].any() and mask[:, 2 * N:].any())
    # isolated pixels (4-neighbourhood, periodic)
    P = np.pad(A, 1, mode="wrap")
    nb4 = P[:-2, 1:-1] + P[2:, 1:-1] + P[1:-1, :-2] + P[1:-1, 2:]
    rep["isolated_pixels"] = int(np.sum((A == 1) & (nb4 == 0)))
    # diagonal-only contacts: solid pairs touching only through a corner
    d1 = (A[:-1, :-1] == 1) & (A[1:, 1:] == 1) & (A[:-1, 1:] == 0) & (A[1:, :-1] == 0)
    d2 = (A[:-1, 1:] == 1) & (A[1:, :-1] == 1) & (A[:-1, :-1] == 0) & (A[1:, 1:] == 0)
    rep["diagonal_only_contacts"] = int(d1.sum() + d2.sum())
    # thinnest-member estimate: 2 x min(distance transform) on the medial ridge
    edt = ndimage.distance_transform_edt(np.tile(A, (3, 3)))[N:2 * N, N:2 * N]
    ridge = np.ones_like(A, dtype=bool)
    Pe = np.pad(ndimage.distance_transform_edt(np.tile(A, (3, 3))), 0)[0:3*N,0:3*N]
    E = Pe[N - 1:2 * N + 1, N - 1:2 * N + 1]
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            if dy == dx == 0:
                continue
            ridge &= E[1 + dy:N + 1 + dy, 1 + dx:N + 1 + dx] <= E[1:N + 1, 1:N + 1]
    ridge &= A == 1
    rep["thinnest_member_px"] = float(2 * edt[ridge].min()) if ridge.any() else 0.0
    return rep


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("image")
    ap.add_argument("-N", type=int, required=True, help="output grid size N x N")
    ap.add_argument("-o", "--output", help="output .txt (default: <image>_<N>x<N>.txt)")
    ap.add_argument("--solid", choices=["dark", "light"], default="dark",
                    help="which tone is solid material (default: dark)")
    ap.add_argument("--threshold", type=float, default=0.5)
    ap.add_argument("--rule", choices=["area", "center"], default="area")
    ap.add_argument("--min-frac", type=float, default=0.5,
                    help="area rule: solid-area fraction for pixel=1 (default 0.5)")
    ap.add_argument("--reference-fraction", type=float, default=None,
                    help="analytic solid fraction to compare against")
    ap.add_argument("--rotate", type=int, choices=[0, 90], default=0,
        help="rotate the image 90 degrees counterclockwise before "
             "processing (a rotation, not a flip: chirality is preserved)")
    ap.add_argument("--preview", help="write a PNG preview of the bitmap")
    a = ap.parse_args(argv)

    B, warn_jpeg, grey = load_binary(a.image, a.solid, a.threshold)
    if a.rotate == 90:
        B = np.rot90(B)
    A = resample(B, a.N, a.rule, a.min_frac)
    out = a.output or f"{os.path.splitext(a.image)[0]}_{a.N}x{a.N}.txt"
    np.savetxt(out, A, fmt="%d")
    if a.preview:
        Image.fromarray(((1 - A) * 255).astype(np.uint8)).save(a.preview)

    rep = checks(A)
    ok = True
    print(f"img2homog report — {a.image} -> {out}  ({a.N}x{a.N}, rule={a.rule})")
    if warn_jpeg:
        print("  [WARN] JPEG input: compression noise near edges may create "
              "spurious pixels; prefer PNG.")
    if grey > 0.05:
        print(f"  [WARN] {grey*100:.1f}% of source pixels are grey (anti-aliasing?); "
              f"binarization threshold = {a.threshold}.")
    print(f"  solid fraction          : {rep['solid_fraction']:.4f}"
          + (f"  (reference {a.reference_fraction:.4f}, "
             f"diff {rep['solid_fraction']-a.reference_fraction:+.4f})"
             if a.reference_fraction is not None else ""))
    if a.reference_fraction is not None and \
       abs(rep['solid_fraction'] - a.reference_fraction) > 0.02:
        print("  [WARN] solid fraction deviates >0.02 from reference.")
    print(f"  single solid component  : {rep['single_component']}"
          f"  (components in cell: {rep['n_components_in_cell']})")
    print(f"  periodic span x / y     : {rep['periodic_span_x']} / {rep['periodic_span_y']}")
    print(f"  isolated pixels         : {rep['isolated_pixels']}")
    print(f"  diagonal-only contacts  : {rep['diagonal_only_contacts']}")
    print(f"  thinnest member (est.)  : {rep['thinnest_member_px']:.1f} px")
    if not (rep['single_component'] and rep['periodic_span_x'] and rep['periodic_span_y']):
        print("  [FAIL] the solid phase is not a single periodically connected "
              "component — the homogenized stiffness will be spurious.")
        ok = False
    if rep['isolated_pixels'] or rep['diagonal_only_contacts']:
        print("  [WARN] isolated pixels / corner-only contacts present: these act "
              "as spurious voids or point hinges (cf. rotating-unit pathology).")
    if 0 < rep['thinnest_member_px'] < 4:
        print("  [WARN] thinnest member < 4 px: increase N (mesh-convergence "
              "criterion: >= 4 elements across the thinnest member).")
    print("  orientation: row 1 = TOP of the image; flipping rows reverses "
          "chirality.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
