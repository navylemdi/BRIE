import yaml
import sys

class Material:
    """Class material
    """
    def __init__(self, name: str):
        """Initialization of the Material class

        Args:
            name (str): Name of the material according to the material database
        """
        self.name = name
        self.data = self.load_data()
        self.E = float(self.data[name]["E"])
        self.Rho = self.data[name]["Rho"]
        self.Nu = self.data[name]["Nu"]
        self.CTE = self.data[name]["CTE"]
        self.HRC_min = self.data[name]["HRC_min"]
        self.HRC_max = self.data[name]["HRC_max"]
        self.T_lim = self.data[name]["T_lim"]

    def load_data(self):
        """A loader for the material data

        Returns:
            dict: Dictionary of the material properties of the selected material
        """
        try:
            with open("Database_material.yaml", 'r') as data:
               temp=yaml.load(data, Loader=yaml.FullLoader)
               temp2=temp[self.name]
               return temp
        except KeyError:
            print('Unknown material')
            print('Try between : \n', 'AISI 440C (M&I)\n',
                                      'AISI 440C (S&T)\n',
                                      'Ceramic\n',
                                      'AMS5898\n',
                                      'AISI M50\n',
                                      'SAE 52100 (M&I)\n',
                                      'SAE 52100 (S&T)\n')
            sys.exit()

    def display_properties(self):
        """Print the material properties of the selected material
        """
        print('Material properties of', self.name,':', self.data[self.name])