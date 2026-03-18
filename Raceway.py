import numpy as np
import Material

class Raceway():
    """Raceway object class
    """
    def __init__(self, Material:Material, type: str, r:float, d:float, b: float, ds : float):
        """Initialization of the Raceway object

        Args:
            Material (Material): Material of the raceway
            type (str): Type of the racweay (Inner or Outer)
            r (float): Raceway groove curvature radius [m]
            d (float): Raceway inner or outer diameter [m]
            b (float): Raceway width [m]
            ds (float): Shoulder diameter [m]
        """
        self.Mat = Material
        self.type = type
        self.check_type()
        self.r = r
        self.d = d
        self.b = b
        self.ds = ds

    def check_type(self):
        if self.type!='Inner' and self.type!='Outer':
            raise Exception("Body type not a Inner or Outer ring, please chose between Inner or Outer")
        else:
            print('Body type correct')
