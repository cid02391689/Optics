import numpy as np
from raytracer.physics import refract

class OpticalElement:       
    def intercept(self, ray):
        """
        This method takes a Ray object as an argument and calculates where its path will
        intercept with the element. If no intercept is found it should return None
        """
        raise NotImplementedError('intercept() needs to be implemented in derived classes')
        
    def propagate_ray(self, ray):
        """
        This method takes a Ray object as an argument and will propagate it through the element. To do this, the method will need to:

        1. Work out if and where the ray's path will intercept with the object.
        2. Calculate any change that the element may cause to the ray's direcion vector.
        3. Update the ray's internal state to reflect this new intercept position and direction.
        """
        raise NotImplementedError('propagate_ray() needs to be implemented in derived classes')
        
class GeneralSurface(OpticalElement):
    """
    A class to represent a general surface optical element.
    """
    
    def __init__(self, z_0, aperture, curvature):
        """
    A general optical surface, either planar or spherical.

    Parameters:
    z_0 : float
        The z-coordinate of the surface vertex.
    aperture : float
        The semi-diameter (clear aperture) of the surface.
    curvature : float
        The curvature (1/R).  Zero => plane, nonzero => sphere of radius R = 1/curvature.

    Attributes:
    _z0 : float
        Stored vertex position along z-axis.
    _aperture : float
        Stored semi-diameter.
    _curvature : float
        Stored curvature value.
    _radius : float
        Sphere radius (∞ for plane).
    """
        if z_0 is None:
            z_0 = 0.0
        if curvature is None:
            curvature = 0.0
        if aperture is None:
            aperture = np.inf
        self._z0 = float(z_0)
        self._aperture = float(aperture)
        self._radius = np.inf if curvature == 0.0 else 1.0 / curvature
        self._curvature = curvature

    def z_0(self):
        """
        Get the z-coordinate of the plane.
        Returns: float
            The z-position of the surface vertex.
        """
        return self._z0
    def aperture(self):
        """
        Returns:float
            The semi-diameter of the surface.

        Get the aperture of the plane.
        """
        return self._aperture
    def curvature(self):
        """
        Get the curvature of the plane.
        Returns: float
            The curvature (1/R). Zero for a plane.
        """
        return self._curvature
    
    def intercept(self, ray):
        """
        Compute the intersection point of a ray with this surface.

        Parameters:
        ray : Ray
            The incident ray.

        Returns:
        np.ndarray or None
            The 3D point of intersection, or None if no valid intercept
            (e.g., misses aperture or is behind the ray origin).
        """
        P = np.array(ray.pos())
        k = np.array(ray.direc())

        # Plane case
        if self._curvature == 0.0:
            z0 = self.z_0()
            l = (z0 - P[2]) / k[2]
            if l <= 0:
                return None
            Q = P + l * k
            return Q

        else:
            # Sphere case: solve the quadratic for ray–sphere intersection            
            sign = 1.0 if self._curvature >= 0 else -1.0
            R = 1.0 / abs(self._curvature)
            O = np.array([0.0, 0.0, self.z_0() + sign * R])
            r = P - O

            if (np.dot(r, k)**2 - (np.dot(r, r) - R**2)) < 0:
                return None 
            l1 = -np.dot(r, k) + np.sqrt((np.dot(r, k))**2 - (np.dot(r, r) - R**2))
            l2 = -np.dot(r, k) - np.sqrt((np.dot(r, k))**2 - (np.dot(r, r) - R**2))
            
            eps = 1e-9
            roots = [l for l in (l1, l2) if l > eps]
            if not roots:
                return None
            if len(roots) == 1:
                l_hit = roots[0]
            else: 
                if self._curvature < 0:
                    l_hit = max(roots)
                else: 
                    l_hit = min(roots)

            Q = P + l_hit * k

            if np.hypot(Q[0]-O[0], Q[1]-O[1]) > self.aperture():
                return None

            return Q

class SphericalRefraction(GeneralSurface):
        """
    A refractive spherical surface that applies Snell's law.

    Parameters:
    z_0 : float
        Vertex z-coordinate of the sphere.
    aperture : float
        Semi-diameter of the clear aperture.
    curvature : float
        Surface curvature (1/R).
    n_1 : float
        Refractive index before the surface.
    n_2 : float
        Refractive index after the surface.

    Attributes:
    __n_1 : float
        Index before refraction.
    __n_2 : float
        Index after refraction.
    """
        
        def __init__(self, z_0, aperture, curvature, n_1, n_2):
            super().__init__(z_0=z_0, aperture=aperture, curvature=curvature)
            self.__n_1 = float(n_1)
            self.__n_2 = float(n_2)


        def n_1(self):
            """
            Get the refractive index of the medium before the element.
            
            Returns:
                float: The refractive index of the medium before the element.
            """
            return self.__n_1

        def n_2(self):
            """
            Get the refractive index of the medium after the element.
            
            Returns:
                float: The refractive index of the medium after the element.
            """
            return self.__n_2

        def propagate_ray(self, ray):
            """
            Propagate the ray through the spherical element.
            Parameters:
            ray : Ray
                The incident ray to be refracted.

            Returns:
                Ray or None
            The updated ray with appended vertex and direction,
            or None if total internal reflection occurs or misses aperture.
            """
            Q = self.intercept(ray)
            if Q is None:
                return None
            
            if self._curvature == 0.0:
                current = ray.direc()
                normal = np.array([0.0, 0.0, -np.sign(current[2])])

            else:
                R = 1.0 / abs(self._curvature)
                sign = 1.0 if self._curvature >= 0 else -1.0
                O = np.array([0.0, 0.0, self.z_0() + sign * R])
                normal = (Q - O) / np.linalg.norm(Q - O)
        
                current = ray.direc()
                if np.dot(current, normal) > 0: #  if the ray is going out of the sphere
                    normal = -normal

            k2 = refract(current, normal, self.n_1(), self.n_2()) # new direction
            if k2 is None:
                return None

            ray.append(Q, k2)
            return ray
        
        def focal_point(self):
            """
            Calculate the focal point of the spherical element.
            
            Returns:
                np.ndarray: The focal point of the spherical element.
            """
            R = 1.0 / self._curvature
            f = self.__n_2 * R / (self.__n_2 - self.__n_1)
            return self.z_0() + f
           
class OutputPlane(GeneralSurface):
    """
        A class to represent an output plane.
        
        Attributes:
            z_0 (float): The z-coordinate of the plane.
    """
        
    def __init__(self, z_0):
        super().__init__(z_0=z_0, aperture=np.inf, curvature=0.0)

    def propagate_ray(self, ray):
        """
        Propagate the ray through the output plane.
        Parameters"
        ray : Ray
            The ray to intercept.

        Returns:
        Ray or None
            The ray with an appended vertex at this plane.
        """
        Q = self.intercept(ray)
        if Q is None:
            return None
        ray.append(pos=Q, direc=ray.direc())
        return ray
        
class SphericalReflection(GeneralSurface):
    """
    A spherical mirror that reflects rays by the law of reflection.

    Parameters:
    z_0 : float
        Vertex z-coordinate of the mirror.
    aperture : float
        Semi-diameter of the mirror.
    curvature : float
        Mirror curvature (1/R). Negative => concave, positive => convex.

    Attributes:
    _radius : float
        Radius of curvature (signed).
    """
    def __init__(self, z_0, aperture, curvature):
        super().__init__(z_0=z_0, aperture=aperture, curvature=curvature)
        self.__n = None

    def propagate_ray(self, ray):
        """
        Propagate a ray by reflecting it off the spherical surface.
        Computes intersection, surface normal, applies reflection formula,
        and appends the reflected segment to the ray.

        Parameters:
        ray : Ray
            The incident ray to reflect.

        Returns:
        Ray or None
            The updated ray with a reflected vertex and direction,
            or None if it misses the mirror aperture.
        """
        Q = self.intercept(ray)
        if Q is None:
            return None
        
        # curvature != 0 
        R = 1.0 / self._curvature
        O = np.array([0.0, 0.0, self.z_0() + R])
        geom = (Q - O)
        normal = geom/np.linalg.norm(geom)
        incident = ray.direc()
        if np.dot(incident, normal) > 0:
            normal = -normal

        # reflection
        rout = incident - 2 * np.dot(incident, normal) * normal
        ray.append(Q, rout)
        return ray

    def focal_point(self):
        """
        Return focal point location along z-axis.
        For a spherical mirror, f = |R|/2, and its sign depends on the curvature:
          - Concave (R<0): real focus at z0 - f
          - Convex (R>0): virtual focus at z0 + f
        """
        f = abs(self._radius) / 2.0

        if self._radius < 0:
            return self._z0 - f
        return self._z0 + f