import numpy as np
from .general import GeneralModel
from src.radiative_shift.atomspecies import AtomSpecies


class DisorderedModel(GeneralModel):

    def __init__(self, length, radius, density, medium_atom: AtomSpecies, reference_atom: AtomSpecies):
        super().__init__()
        self.medium_atom = medium_atom
        self.reference_atom = reference_atom

        n = int(length * (density * np.pi * radius ** 2))

        # Generate random points using numpy's random functions
        phi = np.random.uniform(0, 2 * np.pi, size=n)
        u = np.random.uniform(0, radius ** 2, size=n)
        x = np.sqrt(u) * np.cos(phi)
        y = np.sqrt(u) * np.sin(phi)
        z = np.random.uniform(0, length, size=n)

        self.x = np.array(x)
        self.y = np.array(y)
        self.z = np.array(z)

        self._refresh_properties()

    def _refresh_properties(self):
        self.properties.width = np.amax(np.sqrt([x ** 2 + y ** 2 for x, y in zip(self.x, self.y)]))
        self.properties.length = np.amax(self.z)
        self.properties.n_atoms = len(self.x)
        self.properties.density = len(self.x) / self.properties.length / (np.pi * self.properties.width ** 2)

        self.properties.width /= self.medium_atom.lbar
        self.properties.length /= self.medium_atom.lbar
        self.properties.density *= self.medium_atom.lbar ** 3