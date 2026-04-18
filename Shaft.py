from Material import Material

class Shaft():
    def __init__(self, d: float, tolmin: float, tolmax: float, Material : Material):
        self.d = d
        self.tolmin = tolmin
        self.tolmax = tolmax
        