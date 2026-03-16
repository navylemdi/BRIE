import numpy as np
import Material

class Ball():
    """Ball object class
    """
    def __init__(self, Material: Material, D: float):
        """Initialization of the Ball object

        Args:
            Material (Material): Ball material
            D (float): Ball diameter in meters
        """
        self.Mat = Material
        self.D = D
    
    def mass(self):
        """To calculate the ball's mass

        Returns:
            float: Ball's mass
        """
        return self.Mat.Rho * np.pi * (self.d**3)/6

