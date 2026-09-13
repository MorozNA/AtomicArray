import numpy as np
from .general import GeneralModel
from src.radiative_shift.atomspecies import AtomSpecies


class DisorderedSphere(GeneralModel):

    def __init__(self, density, radius, medium_atom: AtomSpecies, reference_atom: AtomSpecies):
        super().__init__()
        self.medium_atom = medium_atom
        self.reference_atom = reference_atom

        n = int(density * 4 / 3 * np.pi * radius ** 3)
        u = np.random.uniform(-1, 1, size=n)
        phi = np.random.uniform(0, 2 * np.pi, size=n)
        r = (np.random.uniform(0, 1, size=n) ** (1 / 3)) * radius
        self.x = r * np.cos(phi) * (1 - u ** 2) ** 0.5
        self.y = r * np.sin(phi) * (1 - u ** 2) ** 0.5
        self.z = r * u

        self._refresh_properties()
        self.write_log()

    def _refresh_properties(self):
        self.properties.width = np.amax(np.sqrt([x ** 2 + y ** 2 for x, y in zip(self.x, self.y)]))
        self.properties.length = self.properties.width
        self.properties.n_atoms = len(self.x)
        self.properties.density = self.properties.n_atoms / (4 / 3 * np.pi * self.properties.width ** 3)

        self.properties.width /= self.medium_atom.lbar
        self.properties.length /= self.medium_atom.lbar
        self.properties.density *= self.medium_atom.lbar ** 3
