# Optics — 3D Geometric Ray Tracer (Spherical Optics)

A lightweight 3D **geometric optics** ray tracer for **planar / spherical** optical elements, built for **quantitative lens performance analysis** (spot diagrams, RMS spot size) and simple lens design comparisons.

---

## Key features

- **Generalized optical surfaces**
  - Planar and spherical geometry
  - Refraction and reflection
  - Output plane for spot evaluation

- **Physics-correct propagation**
  - Snell’s law refraction
  - Handles **total internal reflection** (TIR) cases
  - Optional consideration of Fresnel loss (kept small for this project)

- **Robust ray–surface intersection**
  - Plane intersection (`z = z0`)
  - Spherical intersection via quadratic solve
  - Physically valid root selection + aperture checking

- **Quantitative evaluation**
  - 3D track plots (ray trajectories)
  - 2D spot diagrams at a chosen plane (typically focal plane)
  - RMS spot size to quantify spherical aberration
  - Comparison against diffraction scale (Δx = λ f / D)

- **Lens experiments included**
  - Plano-convex (PC) vs convex-plano (CP) orientation comparison
  - Thick-lens / principal-plane shift correction (for finite thickness)
  - Simple biconvex optimization exploration

---

## Repository structure

- raytracer (core library: rays, surfaces, optics, analysis)
- tests (pytest)
- README.md
- report
- graph (result from tasks)

## Playbook: tweak parameters to get different behaviours

Most experiments in this repo reduce to:

**(1) ray bundle (input)** + **(2) lens geometry** + **(3) observation plane (where you measure the spot)**

Below is a practical “change X → you get Y” cheat sheet.

### A) Ray bundle (input rays)

Typical parameters you will see in `raytracer/analysis.py` (names may vary):

- `bundle_radius` / `r_max` / `radius`  
  **Increase** → includes more marginal rays → **spherical aberration increases** → RMS spot grows faster.

- `nrings`, `rays_per_ring` (sampling density)  
  **Increase** → smoother spot/RMS curves and less sampling noise, but **slower runtime**.

- `wavelength` / `lambda`  
  Controls the diffraction scale (e.g. \( \Delta x \approx \lambda f / D \)).  
  **Increase** → diffraction-limited spot size increases.

**Try this:** keep the lens fixed and sweep `bundle_radius` from small → large, then plot RMS spot radius.

---

### B) Lens geometry (surfaces + material)

Common lens parameters:

- `R1`, `R2` / `curvature1`, `curvature2`  
  Stronger curvature (smaller \(|R|\) / larger \(|1/R|\)) → **shorter focal length**.  
  Curvature distribution affects aberrations (PC vs CP orientation).

- `thickness`  
  **Increase** → principal plane shift matters more (thick-lens behaviour).

- `aperture` / `D` (clear aperture diameter)  
  **Decrease** → reduces aberration contribution (smaller RMS), but increases diffraction limit (larger \( \Delta x \)).

- `n_lens` (refractive index)  
  **Increase** → stronger refraction → shorter focal length (for the same curvatures).

**Try this:** flip a plano-convex lens orientation (PC ↔ CP) and compare RMS / spot diagrams.

---

### C) Observation plane (where you evaluate the spot)

Typical parameters:

- `z_image` / `z_screen` / output plane position  
  Moving the plane changes the measured spot:  
  sweep around the expected focus → find **minimum RMS** (best focus).

**Try this:** scan `z_image` around the expected focus and plot RMS vs `z_image` to locate best focus numerically.

---

## Quickstart (run locally)

Clone the repo and create a virtual environment:

```bash
git clone https://github.com/cid02391689/Optics.git
cd Optics

python3 -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
pip install numpy matplotlib pytest

