import numpy as np
import Material

class Raceway():
    """Raceway object class
    """
    def __init__(self, Material:Material, r:float, d:float, b: float, ds : float):
        """Initialization of the Raceway object

        Args:
            Material (Material): Material of the raceway
            r (float): Raceway groove curvature radius [m]
            d (float): Raceway inner or outer diameter [m]
            b (float): Raceway width [m]
            ds (float): Shoulder diameter [m]
        """
        self.Mat = Material
        self.r = r
        self.d = d
        self.b = b
        self.ds = ds

    