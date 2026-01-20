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

- raytracer/      # core library (rays, surfaces, optics, analysis)
- tests/          # unit tests (pytest)
- README.md
- report
