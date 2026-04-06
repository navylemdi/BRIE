import numpy as np

class BearingLoads:
    """Classe pour représenter les forces et moments appliqués à un roulement."""

    def __init__(self, Fa: float = 0.0, Fr: float = 0.0, M: float = 0.0):
        """
        Initialise les forces et moments appliqués au roulement.

        Args:
            Fa (float): Force axiale [N].
            Fr (float): Force radiale [N].
            M (float): Moment [Nm].
        """
        self.Fa = Fa  # Force axiale
        self.Fr = Fr  # Force radiale
        self.M = M    # Moment

    def __repr__(self):
        """Représentation textuelle de l'objet."""
        return f"BearingLoads(Fa={self.Fa:.2f} N, Fr={self.Fr:.2f} N, M={self.M:.2f} Nm)"

    def set_loads(self, Fa: float, Fr: float, M: float):
        """Met à jour les forces et moments.

        Args:
            Fa (float): Force axiale [N].
            Fr (float): Force radiale [N].
            M (float): Moment [Nm].
        """
        self.Fa = Fa
        self.Fr = Fr
        self.M = M

    def get_resultant_force(self) -> float:
        """Calcule la force résultante.

        Returns:
            float: Force résultante [N].
        """
        return np.sqrt(self.Fa**2 + self.Fr**2)

    def get_load_angle(self, degrees: bool = True) -> float:
        """Calcule l'angle de la force résultante par rapport à l'axe radial.

        Args:
            degrees (bool): Si True, retourne l'angle en degrés. Sinon, en radians.

        Returns:
            float: Angle de la force résultante.
        """
        angle_rad = np.arctan2(self.Fa, self.Fr)
        return np.degrees(angle_rad) if degrees else angle_rad

    def scale_loads(self, factor: float):
        """Met à l'échelle les forces et moments par un facteur donné.

        Args:
            factor (float): Facteur de mise à l'échelle.
        """
        self.Fa *= factor
        self.Fr *= factor
        self.M *= factor

    def add_loads(self, other: 'BearingLoads'):
        """Ajoute les forces et moments d'un autre objet BearingLoads.

        Args:
            other (BearingLoads): Autre objet BearingLoads.
        """
        self.Fa += other.Fa
        self.Fr += other.Fr
        self.M += other.M

    def to_array(self) -> np.ndarray:
        """Retourne les forces et moments sous forme de tableau NumPy.

        Returns:
            np.ndarray: Tableau contenant [Fa, Fr, M].
        """
        return np.array([self.Fa, self.Fr, self.M])

    def normalize(self, max_load: float = 1.0):
        """Normalise les forces et moments par rapport à une charge maximale.

        Args:
            max_load (float): Charge maximale de référence.
        """
        current_max = max(self.get_resultant_force(), abs(self.M))
        if current_max > 0:
            self.scale_loads(max_load / current_max)

    def is_zero(self, tolerance: float = 1e-6) -> bool:
        """Vérifie si les forces et moments sont nuls (à une tolérance près).

        Args:
            tolerance (float): Tolérance pour considérer une valeur comme nulle.

        Returns:
            bool: True si les forces et moments sont nuls.
        """
        return (abs(self.Fa) < tolerance and
                abs(self.Fr) < tolerance and
                abs(self.M) < tolerance)