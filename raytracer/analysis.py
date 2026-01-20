"""Analysis module."""
import matplotlib.pyplot as plt
from raytracer.rays import Ray, RayBundle
from raytracer.elements import SphericalRefraction, OutputPlane, SphericalReflection
import numpy as np
from raytracer.lenses import PlanoConvex, BiConvex
from scipy.optimize import minimize

def task8():
    """
    Task 8.

    In this function you should check your propagate_ray function properly
    finds the correct intercept and correctly refracts a ray. Don't forget
    to check that the correct values are appended to your Ray object.
    """
    sr = SphericalRefraction(z_0=10, curvature=0.02, n_1=1.0, n_2=1.5, aperture=50.0)
    results = []
    # Propagate rays at different initial x-positions (on-axis and off-axis)
    for x in [-10.0, -5.0, 0.0, 5.0, 10.0]:
        ray = Ray(pos=[x, 0.0, 0.0], direc=[0.0, 0.0, 1.0])
        sr.propagate_ray(ray)
        # After propagation, the last vertex is the intercept
        intercept_point = ray.vertices()[-1]
        # Direction property returns the current direction
        refracted_direction = ray.direc()
        results.append((x, intercept_point, refracted_direction))
    return results


def task10():
    """
    Task 10.

    In this function you should create Ray objects with the given initial positions.
    These rays should be propagated through the surface, up to the output plane.
    You should then plot the tracks of these rays.
    This function should return the matplotlib figure of the ray paths.

    Returns:
        Figure: the ray path plot.
    """
    sr = SphericalRefraction(z_0=100, curvature=0.03, n_1=1.0, n_2=1.5, aperture=34.0)
    op = OutputPlane(z_0=250)
    fig, ax = plt.subplots(figsize=(8, 6))
    for y in [-20.0, -10.0, -4.0, -1.0, -0.2, 0.0, 0.2, 1.0, 4.0, 10.0, 20.0]:
        ray= Ray(pos=[0.0, y, 0.0], direc=[0.0, 0.0, 1.0])
        sr.propagate_ray(ray)
        op.propagate_ray(ray)
        path = ray.vertices()  # list of 3D points
        zs = [pt[2] for pt in path]
        ys = [pt[1] for pt in path]
        ax.plot(zs, ys, 'o-')

    ax.set_xlabel("z (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title("Ray paths through single spherical surface")
    ax.grid(True)
    return fig


def task11():
    """
    Task 11.

    In this function you should propagate the three given paraxial rays through the system
    to the output plane and the tracks of these rays should then be plotted.
    This function should return the following items as a tuple in the following order:
    1. the matplotlib figure object for ray paths
    2. the calculated focal point.

    Returns:
        tuple[Figure, float]: the ray path plot and the focal point
    """
    sr = SphericalRefraction(z_0=100, curvature=0.03, n_1=1.0, n_2=1.5, aperture=34.0)
    fz = sr.focal_point()
    op = OutputPlane(z_0=fz)
    fig = plt.figure(figsize=(8,6))
    ax  = fig.add_subplot(111, projection='3d')
    for pos in [[0.1, 0.1, 0.0], [0.0, 0.0, 0.0], [-0.1, -0.1, 0.0]]:
        ray = Ray(pos=pos, direc=[0.0, 0.0, 1.0])
        sr.propagate_ray(ray)
        op.propagate_ray(ray)
        path = ray.vertices()
        zs = [pt[2] for pt in path]
        ys = [pt[1] for pt in path]
        xs = [pt[0] for pt in path]
        ax.plot(xs, ys, zs, 'o-')
    ax.set_xlabel("X (mm)")
    ax.set_ylabel("Y (mm)")
    ax.set_title("Ray paths through single spherical surface")
    ax.grid(True)
    return fig , fz

def task12():
    """
    Task 12.

    In this function you should create a RayBunble and propagate it to the output plane
    before plotting the tracks of the rays.
    This function should return the matplotlib figure of the track plot.

    Returns:
        Figure: the track plot.
    """
    sr = SphericalRefraction(z_0=100, curvature=0.03, n_1=1.0, n_2=1.5, aperture=34.0)
    fz = sr.focal_point()
    op = OutputPlane(z_0=fz)
    rb = RayBundle()
    rb.propagate_bundle([sr, op])
    fig = rb.track_plot()
    return fig


def task13():
    """
    Task 13.

    In this function you should again create and propagate a RayBundle to the output plane
    before plotting the spot plot.
    This function should return the following items as a tuple in the following order:
    1. the matplotlib figure object for the spot plot
    2. the simulation RMS

    Returns:
        tuple[Figure, float]: the spot plot and rms
    """
    sr = SphericalRefraction(z_0=100, curvature=0.03, n_1=1.0, n_2=1.5, aperture=34.0)
    fz = sr.focal_point()
    op = OutputPlane(z_0=fz)
    rb = RayBundle()
    rb.propagate_bundle([sr, op])
    fig = rb.spot_plot()
    rms = rb.rms()
    return fig, rms

def task14():
    """
    Task 14.

    In this function you will trace a number of RayBundles through the optical system and
    plot the RMS and diffraction scale dependence on input beam radii.
    This function should return the following items as a tuple in the following order:
    1. the matplotlib figure object for the diffraction scale plot
    2. the simulation RMS for input beam radius 2.5
    3. the diffraction scale for input beam radius 2.5

    Returns:
        tuple[Figure, float, float]: the plot, the simulation RMS value, the diffraction scale.
    """
    sr = SphericalRefraction(z_0=100, curvature=0.03, n_1=1.0, n_2=1.5, aperture=34.0)
    fz = sr.focal_point()
    f = fz - sr.z_0()
    op = OutputPlane(z_0=fz)
    wavelength = 588e-6
    radii = np.linspace(0.1, 10.0, 100)
    rms_list = []
    diffraction_scale_list = []
    for r in radii:
        rb = RayBundle(rmax=r, nrings=5, multi=6)
        rb.propagate_bundle([sr, op])
        rms = rb.rms()
        rms_list.append(rms)
        diffraction_scale = wavelength * f / (2 * r)
        diffraction_scale_list.append(diffraction_scale)
    
    idx25 = np.argmin(np.abs(radii - 2.5)) # Find index of closest value to 2.5
    rms_25 = rms_list[idx25]
    dx_25  = diffraction_scale_list[idx25]
    # Plotting
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(radii, rms_list, label="Geometric RMS")
    ax.plot(radii, diffraction_scale_list,  '--', label="Diffraction Δx")
    ax.set_xlabel("Beam radius R (mm)")
    ax.set_ylabel("RMS/Diffraction scale (mm)")
    ax.set_title("RMS and Diffraction scale vs Beam radius")
    ax.legend()
    ax.grid(True)    
    return fig, rms_25, dx_25


def task15():
    """
    Task 15.

    In this function you will create plano-convex lenses in each orientation and propagate a RayBundle
    through each to their respective focal point. You should then plot the spot plot for each orientation.
    This function should return the following items as a tuple in the following order:
    1. the matplotlib figure object for the spot plot for the plano-convex system
    2. the focal point for the plano-convex lens
    3. the matplotlib figure object for the spot plot for the convex-plano system
    4  the focal point for the convex-plano lens


    Returns:
        tuple[Figure, float, Figure, float]: the spot plots and rms for plano-convex and convex-plano.
    """
    # Create plano-convex lens with curved side first
    pc = PlanoConvex(z_0=100, curvature=None, curvature1=0., curvature2=-0.02, n_inside=1.5168, n_outside=1., thickness=5., aperture=50.)
    fz_pc = pc.focal_point()
    outputplane_pc = OutputPlane(z_0=fz_pc)
    rb_pc = RayBundle()
    rb_pc.propagate_bundle([pc._surface1, pc._surface2, outputplane_pc])
    fig_pc = rb_pc.spot_plot()
    ax_pc = fig_pc.axes[0]                    
    ax_pc.set_title("Plano-Convex Lens Spot Plot")
    ax_pc.set_xlabel("x (mm)")
    ax_pc.set_ylabel("y (mm)")
    ax_pc.grid(True)

    # Create convex-plano lens with curved side second
    cp = PlanoConvex(z_0=100, curvature=None, curvature1=0.02, curvature2=0., n_inside=1.5168, n_outside=1., thickness=5., aperture=50.)
    fz_cp = cp.focal_point()
    outputplane_cp = OutputPlane(z_0=fz_cp)
    rb_cp = RayBundle()
    rb_cp.propagate_bundle([cp._surface1, cp._surface2, outputplane_cp])
    fig_cp = rb_cp.spot_plot()
    ax_cp = fig_cp.axes[0]
    ax_cp.set_title("Convex-Plano Lens Spot Plot")
    ax_cp.set_xlabel("x (mm)")
    ax_cp.set_ylabel("y (mm)")
    ax_cp.grid(True)
   
    return fig_pc, fz_pc, fig_cp, fz_cp


def task16():
    """
    Task 16.

    In this function you will be again plotting the radial dependence of the RMS and diffraction values
    for each orientation of your lens.
    This function should return the following items as a tuple in the following order:
    1. the matplotlib figure object for the diffraction scale plot
    2. the RMS for input beam radius 3.5 for the plano-convex system
    3. the RMS for input beam radius 3.5 for the convex-plano system
    4  the diffraction scale for input beam radius 3.5

    Returns:
        tuple[Figure, float, float, float]: the plot, RMS for plano-convex, RMS for convex-plano, diffraction scale.
    """
    wavelength = 588e-6
    radii = np.linspace(0.1, 10.0, 100)
    rms_pc_list = []
    rms_cp_list = []
    diff_list = []
    for r in radii:
        # Plano-Convex lens
        pc = PlanoConvex(z_0=100, curvature=None, curvature1=0., curvature2=-0.02, n_inside=1.5168, n_outside=1., thickness=5., aperture=50.)
        fz_pc = pc.focal_point()
        outputplane_pc = OutputPlane(z_0=fz_pc)
        rb_pc = RayBundle(rmax=r, nrings=5, multi=6)
        rb_pc.propagate_bundle([pc._surface1, pc._surface2, outputplane_pc])
        rms_pc = rb_pc.rms()
        rms_pc_list.append(rms_pc)

        # Convex-Plano lens
        cp = PlanoConvex(z_0=100, curvature=None, curvature1=0.02, curvature2=0., n_inside=1.5168, n_outside=1., thickness=5., aperture=50.)
        fz_cp = cp.focal_point()
        outputplane_cp = OutputPlane(z_0=fz_cp)
        rb_cp = RayBundle(rmax=r, nrings=5, multi=6)
        rb_cp.propagate_bundle([cp._surface1, cp._surface2, outputplane_cp])
        rms_cp = rb_cp.rms()
        rms_cp_list.append(rms_cp)

        # diffraction scale
        f = pc.focal_point() - pc.z_0() - pc.thickness()
        diffraction_scale =  wavelength * f / (2 * r)
        diff_list.append(diffraction_scale)


    # plotting
    fig, ax = plt.subplots()
    ax.plot(radii, rms_pc_list, label="Plano→Convex RMS")
    ax.plot(radii, rms_cp_list, label="Convex→Plano RMS")
    ax.plot(radii, diff_list, label = "diffraction scale")
    ax.set_xlabel("Beam radius w (mm)")
    ax.set_ylabel("Spot RMS (mm)")
    ax.legend()
    ax.grid(True, which="both", ls=":")

    # Find index of closest value to 3.5
    r0 = 3.5
    idx = np.abs(radii - r0).argmin()
    rms_pc_3p5 = rms_pc_list[idx]
    rms_cp_3p5 = rms_cp_list[idx]
    diff_3p5 = diff_list[idx]

    return fig, rms_pc_3p5, rms_cp_3p5, diff_3p5 


def task17():
    """
    Task 17.

    In this function you will be first plotting the spot plot for your PlanoConvex lens with the curved
    side first (at the focal point). You will then be optimising the curvatures of a BiConvex lens
    in order to minimise the RMS spot size at the same focal point. This function should return
    the following items as a tuple in the following order:
    1. The comparison spot plot for both PlanoConvex (curved side first) and BiConvex lenses at PlanoConvex focal point.
    2. The RMS spot size for the PlanoConvex lens at focal point
    3. the RMS spot size for the BiConvex lens at PlanoConvex focal point

    Returns:
        tuple[Figure, float, float]: The combined spot plot, RMS for the PC lens, RMS for the BiConvex lens
    """
    # curved-first planoconvex lens
    cpl = PlanoConvex(z_0=100, curvature=None, curvature1=0.02, curvature2=0., n_inside=1.5168, n_outside=1., thickness=5., aperture=50.)
    fz_cpl = cpl.focal_point()
    outputplane_cp = OutputPlane(z_0=fz_cpl)
    rb_cpl = RayBundle()
    rb_cpl.propagate_bundle([cpl._surface1, cpl._surface2, outputplane_cp])
    rms_cpl = rb_cpl.rms()

    # biconvex
    # Define necessary parameters for the BiConvex lens and ray bundle
    z0 = 100
    n_inside = 1.5168
    n_outside = 1.0
    thickness = 5.0
    aperture = 50.0
    beam_radius = 5.0
    z_focus = fz_cpl  # Use the focal point of the plano-convex lens for comparison

    def objective(curv):
        c1, c2 = curv
        bi = BiConvex(z_0=z0, curvature1=c1, curvature2=c2, n_inside=n_inside, n_outside=n_outside, thickness=thickness, aperture=aperture)
        rb = RayBundle()
        rb.propagate_bundle(bi.surfaces() + [OutputPlane(z_0=z_focus)])
        return rb.rms()
    
    initial_guess = np.array([+0.02, -0.02])
    bounds = [(-0.1, 0.1), (-0.1, 0.1)]   
    result = minimize(fun=objective,x0=initial_guess,bounds=bounds,method='L-BFGS-B',options={'ftol':1e-6, 'disp': True})

    c1_opt, c2_opt = result.x
    rms_bi = result.fun
    bi = BiConvex(z_0=100., curvature1=c1_opt, curvature2=c2_opt,n_inside=1.5168, n_outside=1.,thickness=5., aperture=50.)
    rb_bi = RayBundle()
    rb_bi.propagate_bundle(bi.surfaces() + [OutputPlane(z_0=z_focus)])

    fig, ax = plt.subplots(figsize=(8,6))
    ax.set_aspect('equal')

    # PlanoConvex
    for ray in rb_cpl:
        x, y, _ = ray.vertices()[-1]
        ax.plot(x, y, 'o', color='C0', alpha=0.6, label='_nolegend_')
    # BiConvex
    for ray in rb_bi:
        x, y, _ = ray.vertices()[-1]
        ax.plot(x, y, 'o', color='C1', alpha=0.6, label='_nolegend_')

    ax.plot([], [], 'o', color='C0', alpha=0.6, label='PlanoConvex')
    ax.plot([], [], 'o', color='C1', alpha=0.6, label='BiConvex')
    ax.set_xlabel("X (mm)")
    ax.set_ylabel("Y (mm)")
    ax.set_title(f"Spot Plot Comparison at z={z_focus:.2f} mm")
    ax.legend(loc='best')
    ax.grid(True)

    rms_bi = rb_bi.rms()
    return fig, rms_cpl, rms_bi

def task18():
    """
    Task 18.

    In this function you will be testing your reflection modelling. Create a new SphericalReflecting surface
    and trace a RayBundle through it to the OutputPlane.This function should return
    the following items as a tuple in the following order:
    1. The track plot showing reflecting ray bundle off SphericalReflection surface.
    2. The focal point of the SphericalReflection surface.

    Returns:
        tuple[Figure, float]: The track plot and the focal point.

    """
    sph_refl = SphericalReflection(z_0=100., aperture=6., curvature=-0.02)
    op = OutputPlane(z_0=50.0)
    rb = RayBundle()
    rb.propagate_bundle([sph_refl, op])
    
    fig = rb.track_plot()  
    fig.axes[0].set_title("Ray Trace for Spherical Reflection")
    fig.axes[0].set_xlabel("x (mm)")
    fig.axes[0].set_ylabel("y (mm)")

    z_focus = sph_refl.focal_point()

    return fig, z_focus


if __name__ == "__main__":

    # Run task 8 function
    task8()

    # Run task 10 function
    # FIG10 = task10()

    # Run task 11 function
    # FIG11, FOCAL_POINT = task11()

    # Run task 12 function
    # FIG12 = task12()

    # Run task 13 function
    # FIG13, TASK13_RMS = task13()

    # Run task 14 function
    # FIG14, TASK14_RMS, TASK14_DIFF_SCALE = task14()

    # Run task 15 function
    # FIG15_PC, FOCAL_POINT_PC, FIG15_CP, FOCAL_POINT_CP = task15()

    # Run task 16 function
    # FIG16, PC_RMS, CP_RMS, TASK16_DIFF_SCALE = task16()

    # Run task 17 function
    # FIG17, CP_RMS, BICONVEX_RMS = task17()

    # Run task 18 function
    # FIG18, FOCAL_POINT = task18()

    plt.show()
