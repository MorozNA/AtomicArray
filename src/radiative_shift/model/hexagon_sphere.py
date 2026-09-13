import numpy as np
from .hexagon import HexagonModel
from src.radiative_shift.atomspecies import AtomSpecies


class HexagonSphere(HexagonModel):

    def __init__(self, radius, density, medium_atom: AtomSpecies, reference_atom: AtomSpecies):
        super().__init__(2 * radius, radius, density, medium_atom, reference_atom)
        excess = []
        zm = np.average(self.z)
        for i, (xi, yi, zi) in enumerate(zip(self.x, self.y, self.z)):
            if xi ** 2 + yi ** 2 + (zi - zm) ** 2 > radius ** 2:
                excess.append(i)

        self.x = np.delete(self.x, excess)
        self.y = np.delete(self.y, excess)
        self.z = np.delete(self.z, excess)

        self._refresh_properties()
        self.write_log()

    def _refresh_properties(self):
        self.properties.width = np.amax(np.sqrt([x ** 2 + y ** 2 for x, y in zip(self.x, self.y)]))
        self.properties.length = self.properties.width
        self.properties.n_atoms = len(self.x)
        self.properties.density = len(self.x) / (4 / 3 * np.pi * self.properties.width ** 3)

        self.properties.width /= self.medium_atom.lbar
        self.properties.length /= self.medium_atom.lbar
        self.properties.density *= self.medium_atom.lbar ** 3
