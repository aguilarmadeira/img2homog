# img2homog

Convert a black-and-white image of a periodic unit cell into an `N x N`
0/1 bitmap ready for finite-element **asymptotic homogenization**
(Guedes & Kikuchi, 1990), with the validation checks that decide whether
the bitmap is actually usable.

Converting an image to 0/1 is trivial; converting it *well* is not. A
naive conversion silently produces bitmaps whose homogenized stiffness is
meaningless: disconnected solid phases, single-pixel point hinges,
anti-aliasing noise, walls thinner than the mesh can resolve, or a
flipped chiral cell with the wrong handedness. `img2homog` performs the
conversion **and prints a validation report** covering exactly these
failure modes.


## Web version (no installation)

A browser version is available at
**<https://aguilarmadeira.github.io/img2homog/>**: drop the image, set the
options, read the validation report and download the `.txt` matrix. It runs
the repository's own Python in your browser (via Pyodide), so results are
identical to the command line and the image never leaves your computer.
The first visit takes 10-20 s to load the Python environment.

## Usage

```bash
python img2homog.py cell.png -N 200 -o cell_200.txt
python img2homog.py cell.png -N 100 --solid dark --threshold 0.5 \
    --reference-fraction 0.345 --preview cell_100.png
```

Output: an `N x N` matrix of space-separated 0/1 values (1 = solid,
0 = void), one matrix row per line — readable in MATLAB with
`A = load('cell_200.txt');`.

### Options

| option | meaning | default |
|---|---|---|
| `-N` | output grid size (required) | — |
| `-o` | output `.txt` path | `<image>_<N>x<N>.txt` |
| `--solid {dark,light}` | which tone is solid material | `dark` |
| `--threshold` | grey-level binarization threshold in [0,1] | `0.5` |
| `--rule {area,center}` | pixel = 1 by solid-area fraction, or by pixel-centre membership | `area` |
| `--min-frac` | area rule: solid fraction for pixel = 1 | `0.5` |
| `--reference-fraction` | analytic solid fraction to compare against | — |
| `--preview` | write a PNG preview of the bitmap | — |

## The validation report

```
img2homog report — examples/missing_rib_v1.png -> missing_rib_v1_100x100.txt  (100x100, rule=area)
  solid fraction          : 0.3468  (reference 0.3450, diff +0.0018)
  single solid component  : True  (components in cell: 1)
  periodic span x / y     : True / True
  isolated pixels         : 0
  diagonal-only contacts  : 0
  thinnest member (est.)  : 12.0 px
```

* **single solid component / periodic span** — the solid phase must be one
  component that connects across the cell boundaries in both directions
  (checked with periodic 4-connectivity on a 3x3 tiling). If not, the
  tool exits with `[FAIL]` and a nonzero status.
* **isolated pixels / diagonal-only contacts** — corner-only contacts act
  as spurious point hinges: they either fuse or disconnect depending on
  the mesh, so the homogenized response does not converge with
  resolution (the classic rotating-squares vertex-hinge pathology).
* **thinnest member** — estimated via the distance-transform medial
  ridge; a warning is issued below 4 px, mirroring the usual
  mesh-convergence criterion of >= 4 elements across the thinnest member.
* **JPEG / anti-aliasing warnings** — JPEG compression and anti-aliased
  edges create grey pixels whose binarization is threshold-sensitive;
  prefer clean PNG input.

## Orientation convention (read this if your cell is chiral)

Row 1 of the output corresponds to the **top** of the image; column 1 to
the left. **Flipping rows reverses chirality**: for a chiral cell this
silently changes the sign of the normal–shear coupling. Keep the
convention consistent between this tool and your homogenization code.

## Examples

`examples/` contains the four planar auxetic cells (re-entrant,
missing-rib, tetrachiral, rotating-unit; wall 0.12 L) used in the
companion study, rendered at 400x400 from their exact analytic geometry.
Their exact solid fractions, for `--reference-fraction`:

| cell | exact solid fraction |
|---|---|
| `reentrant_v1.png` | 0.398 |
| `missing_rib_v1.png` | 0.345 |
| `chiral_v1.png` | 0.289 |
| `rotating_v1.png` | 0.397 |

## Install / test

Requires Python >= 3.9 with `numpy`, `scipy`, `pillow`
(`pip install -r requirements.txt`). Tests: `python tests/test_img2homog.py`.

## Reference

J. M. Guedes, N. Kikuchi, *Preprocessing and postprocessing for materials
based on the homogenization method with adaptive finite element methods*,
Comput. Methods Appl. Mech. Engrg. 83 (1990) 143–198.

## License

MIT — see `LICENSE`.
