"""
Module: rays.py
"""
import numpy as np
from raytracer.genpolar import rtrings
from raytracer.elements import OpticalElement
import matplotlib.pyplot as plt

class Ray:
    """
    Represent a single optical ray in 3D.
    A Ray maintains its current position and direction, and a history of
    all vertex positions visited as it propagates through optical elements.

    Attributes:
    __pos : np.ndarray
        Current 3D position of the ray (shape (3,)).
    __direc : np.ndarray
        Current unit direction vector of the ray (shape (3,)).
    __points : list of np.ndarray
        Ordered list of all positions the ray has visited.
    """

    def __init__(self, pos=None, direc=None):
        """
        Initialize a Ray.

        Parameters :
        pos : 3-D array
            Initial 3D position. Defaults to (0, 0, 0).
        direc : 3-D array
            Initial direction vector. Will be normalized. Defaults to (0, 0, 1).

        Raises:
        TypeError
            If pos or direc cannot be converted to a float array of shape (3,).
        ValueError
            If direction vector has zero length.
        """

        if pos is None:
            self.__pos = np.array([0.0, 0.0, 0.0])
        else:
            try:
                arr = np.array(pos, dtype=float)
            except Exception:
                raise TypeError("pos must be array-like of length 3")
            if arr.shape != (3,):
                raise ValueError("pos must have length 3")
            self.__pos = arr

        if direc is None:
            d = np.array([0.0, 0.0, 1.0])
        else:
            try:
                d = np.array(direc, dtype=float)
            except Exception:
                raise TypeError("direc must be array-like of length 3")
            if d.shape != (3,):
                raise ValueError("direc must have length 3")
        norm = np.linalg.norm(d)
        if norm == 0:
            raise ValueError("Direction vector must be non-zero")
        self.__direc = d / norm
        # Initialize the list of points with the initial position       
        self.__points = [self.__pos.copy()]

    def pos(self):
        """
        Get the position of the ray.

        Returns:
            np.ndarray: The position of the ray.
        """
        return np.array(self.__pos.copy())
    
    def direc(self):
        """
        Get the direction of the ray.

        Returns:
            np.ndarray: The direction of the ray.
        """
        return np.array(self.__direc.copy())
    
    def append(self, pos, direc):
        
        """
        Append a new vertex to the ray, updating position and direction.

        Parameters: 
        pos : array-like of length 3
            New position to append.
        direc : array-like of length 3
            New direction vector (will be normalized).

        Raises: 
        TypeError
            If pos or direc cannot be converted to float arrays.
        ValueError
            If direc has zero length or wrong shape.
        """

        # check append position
        try:
            p = np.array(pos, dtype=float)
        except Exception:
            raise TypeError("pos must be array-like of length 3")
        if p.shape != (3,):
            raise ValueError("pos must have length 3")

        # check append direction
        try:
            d = np.array(direc, dtype=float)
        except Exception:
            raise TypeError("direc must be array-like of length 3")
        if d.shape != (3,):
            raise ValueError("direc must have length 3")
        norm = np.linalg.norm(d)
        if norm == 0:
            raise ValueError("direction vector must be non-zero")

        self.__pos = p
        self.__direc = d / norm
        self.__points.append(self.__pos.copy())

    def vertices(self):
        """
        return: all the points along the ray in the form of a list containing the numpy array position vectors.
        """
        return [pt.copy() for pt in self.__points]
    
class RayBundle:
    """
    A class to represent a bundle of rays in 3D space.

    Attributes:
        __rays (list of Ray): List of Ray objects in the bundle.
    """

    def __init__(self, rmax=5.0, nrings=5, multi=6):
        """
        Initialize the RayBundle by sampling rays in the input plane.

        Parameters:
        rmax : float
            Maximum radius of the outer ring (semi-diameter in mm).
        nrings : int
            Number of concentric rings of rays to sample.
        multi : int
            Number of equally spaced rays per ring.
        """
        self.__rmax = rmax
        self.__nrings = int(nrings)
        self.__multi = multi
        self.__rays = []
        for r, theta in rtrings(rmax=self.__rmax, nrings=self.__nrings,multi=self.__multi):
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            ray = Ray(pos=[x, y, 0.0], direc=[0.0, 0.0, 1.0])
            self.__rays.append(ray)

       
    def propagate_bundle(self, elements):
        """
        Propagate the bundle of rays through optical elements.
        Args:
            element (OpticalElement): A list of optical elements to propagate the rays through.
        """
        for ray in self.__rays:
            for element in elements:               
                element.propagate_ray(ray)

    def track_plot(self):
        """
        Plot the 3D trajectories of all rays in the bundle.

        Returns:
        matplotlib.figure.Figure
            The figure object containing the 3D plot of ray paths.
        """
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection='3d')
        for ray in self.__rays:
            pts = ray.vertices()  
            xs = [p[0] for p in pts]
            ys = [p[1] for p in pts]
            zs = [p[2] for p in pts]
            ax.plot(xs, ys, zs, '-o')
        ax.set_xlabel("X (mm)")
        ax.set_ylabel("Y (mm)")
        ax.set_title("Ray bundle tracks")
        ax.grid(True)
        return fig
    
    def rms(self):
        """
        Calculate the RMS of the ray bundle.
        Returns:
            float: The RMS of the ray bundle.
        """
        x = np.array([ray.pos()[0] for ray in self.__rays])
        y = np.array([ray.pos()[1] for ray in self.__rays])
        rms = np.sqrt(np.mean(x**2 + y**2))
        return rms
    
    def spot_plot(self):
        """
        Create a spot plot of the ray bundle.
        Returns:
            fig: The matplotlib figure object for the spot plot.   
        """
        
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)

        for ray in self.__rays:
            x, y, _ = ray.vertices()[-1]
            ax.plot(x, y, 'o', linestyle='') 
        ax.set_xlabel("X (mm)")
        ax.set_ylabel("Y (mm)")
        ax.set_title("Ray bundle spot plot")
        ax.grid(True)
        return fig   

    @property
    def rays(self):
       return list(self.__rays)

    # for loop can then be applied to raybundle
    def __iter__(self):
        return iter(self.__rays)             