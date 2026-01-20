from raytracer.elements import SphericalRefraction, OutputPlane, OpticalElement

class PlanoConvex(OpticalElement):
    
    def __init__(self, z_0, curvature, curvature1, curvature2, n_inside, n_outside, thickness: float, aperture: float):
        
        """
    One of the optical elements, consisting of one plane surface and one spherical surface.

    Parameters:
    z_0 : float
        z-coordinate of the vertex of the first (left-hand) surface.
    curvature : Optional[float]
        If not None, a single curvature to apply to one side: 
        positive for a concave-toward-source surface, negative for convex.
    curvature1 : float
        Curvature (1/R) of the first surface, if `curvature` is None.
    curvature2 : float
        Curvature (1/R) of the second surface, if `curvature` is None.
    n_inside : float
        Refractive index of the lens material.
    n_outside : float
        Refractive index of the surrounding medium (usually air = 1.0).
    thickness : float
        Center-thickness of the lens (distance between the two vertex planes).
    aperture : float
        Semi-diameter of the clear aperture.

    Attributes:
    _surface1, _surface2 : SphericalRefraction
        The two surface elements composing the lens.
    """

        if curvature is None:
            self._curvature1 = curvature1
            self._curvature2 = curvature2

        else:
            if curvature > 0:  
                self._curvature1 = curvature
                self._curvature2 = 0.0
            else:
                self._curvature1 = 0.0
                self._curvature2 = curvature    

        self.__z_0 = z_0
        self.__n_outside = n_outside
        self.__n_inside = n_inside
        self.__thickness = thickness
        self.__aperture = aperture
        self._z_1 = z_0 + self.__thickness
        
        self._surface1 = SphericalRefraction(z_0=z_0,aperture=aperture,curvature=self._curvature1,n_1=n_outside,n_2=n_inside,)
        self._surface2 = SphericalRefraction(z_0=self._z_1,aperture=aperture,curvature=self._curvature2,n_1=n_inside,n_2=n_outside,)

    def z_0(self):
        """
        Returns: float
        The z-coordinate of the first surface vertex.
        """
        return self.__z_0
    
    def thickness(self):
        """
        Returns:float
        The center thickness of the lens.
        """
        return self.__thickness

    def focal_point(self):
        """
        Calculate the global z-coordinate of the focal point for the lens, using the thick lens formula.
        Returns: float
        The z position of the focal point.
        """
        n, t = self.__n_inside, self.__thickness
        Δn = n - 1.0
        C1, C2 = self._curvature1, self._curvature2
        inv_f = Δn * (C1 - C2 + (Δn * t * C1 * C2) / n)
            
        if inv_f == 0:
            return float('inf')
        f = 1.0 / inv_f

        if C1 == 0.0:
            return self.__z_0 + f + t
        else:
            H = Δn * t * C1 / (n * (C1 - C2))
            return self.__z_0 + f + H
                
    def propagate(self, ray):
        """
        Propagate the ray through the lens. This method applies Snell's law at the first surface, then again
        at the second surface, appending the new segments to `ray`.
        
        Args:
            ray (Ray): The ray to propagate.
        
        Returns:
            Ray: The propagated ray.
        """
        # Propagate through the first surface
        self._surface1.propagate_ray(ray)
        # Propagate through the second surface
        self._surface2.propagate_ray(ray)
        return ray
    
ConvexPlano = PlanoConvex # alias for backward compatibility

class BiConvex(OpticalElement):
    def __init__(self, z_0, curvature1, curvature2, n_inside, n_outside, thickness, aperture):
        
        """
    A bi-convex lens, consisting of two spherical surfaces.

    Parameters:
    z_0 : float
        z-coordinate of the vertex of the first surface.
    curvature1 : float
        Curvature (1/R) of the first (left-hand) spherical surface.
    curvature2 : float
        Curvature (1/R) of the second (right-hand) spherical surface.
    n_inside : float
        Refractive index of the lens material.
    n_outside : float
        Refractive index of the surrounding medium.
    thickness : float
        Center-thickness of the lens (distance between the two vertices).
    aperture : float
        Semi-diameter of the clear aperture.

    Attributes:
    _surface1, _surface2 : SphericalRefraction
        The two spherical surface elements composing the lens.
    """

        self.__z_0 = z_0
        self._curvature1 = curvature1
        self._curvature2 = curvature2
        self.__n_inside = n_inside
        self.__n_outside = n_outside
        self.__thickness = thickness
        self.__aperture = aperture

        self._surface1 = SphericalRefraction(z_0=z_0, curvature=curvature1, n_1=n_outside, n_2=n_inside, aperture=aperture)
        self._surface2 = SphericalRefraction(z_0=z_0 + thickness,curvature=curvature2,n_1=n_inside, n_2=n_outside,aperture=aperture)

    def surfaces(self):
        """
        Returns: list of OpticalElement
            The ordered list of surface elements [first_surface, second_surface]
            that a ray must traverse when passing through this lens.
        """
        return [self._surface1, self._surface2]