import numpy as np
from .general import GeneralModel
from src.radiative_shift.atomspecies import AtomSpecies


class EmptyModel(GeneralModel):

    def __init__(self, medium_atom: AtomSpecies, reference_atom: AtomSpecies):
        super().__init__()
        self.medium_atom = medium_atom
        self.reference_atom = reference_atom

        self._refresh_properties()

    def _refresh_properties(self):
        # TODO: write normal method
        # self.properties.width = np.amax(np.sqrt([x ** 2 + y ** 2 for x, y in zip(self.x, self.y)])) / self.medium_atom.lbar
        # self.properties.length = np.amax(self.z) / self.medium_atom.lbar
        self.properties.width = 1
        self.properties.length = 1
        self.properties.n_atoms = len(self.x)
        # density = n0 * self.medium_atom.lbar ** 3
        self.properties.density = len(self.x) / self.properties.length / (self.properties.width ** 2)
