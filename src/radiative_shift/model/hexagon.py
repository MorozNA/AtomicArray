import numpy as np
from .general import GeneralModel
from src.radiative_shift.atomspecies import AtomSpecies


class HexagonModel(GeneralModel):

    def __init__(self, length, radius, density, medium_atom: AtomSpecies, reference_atom: AtomSpecies):
        super().__init__()
        self.medium_atom = medium_atom
        self.reference_atom = reference_atom

        self.unitradius = density ** (-1 / 3)
        layers = max(1, round(radius / self.unitradius))
        copies = max(1, round(length / self.unitradius))
        n = round(length * (density * np.pi * radius ** 2))
        while True:
            if n >= copies * (1 + 3 * (layers + 1) * (layers + 2)):
                layers += 1
                self.unitradius = radius / layers
                copies = max(1, round(length / self.unitradius))
            else:
                break

        self.unitradius = radius / layers
        copies = max(1, round(length / self.unitradius))
        self._nl = 0

        self.x = [0]
        self.y = [0]
        self.z = [0]

        for _ in range(layers):
            self._add_layer()
        for _ in range(copies):
            self._add_copy()

        self.x = np.array(self.x)
        self.y = np.array(self.y)
        self.z = np.array(self.z)

        self._refresh_properties()
        self.write_log()

    def _add_layer(self):
        self._nl += 1
        layerX = []
        layerY = []
        layerZ = []
        for i in range(6):
            x = self._nl * self.unitradius * np.sin(i * 2 * np.pi / 6)
            y = self._nl * self.unitradius * np.cos(i * 2 * np.pi / 6)
            z = 0
            layerX.append(x)
            layerY.append(y)
            layerZ.append(z)
            for j in range(1, self._nl):
                layerX.append(x + j * self.unitradius * np.sin((i + 2) * 2 * np.pi / 6))
                layerY.append(y + j * self.unitradius * np.cos((i + 2) * 2 * np.pi / 6))
                layerZ.append(z)
        self.x += layerX
        self.y += layerY
        self.z += layerZ

    def _add_copy(self):
        n = round(1 + self._nl * (6 + 6 * self._nl) / 2)
        self.z += n * [self.z[-1] + self.unitradius]
        self.x += self.x[0:n]
        self.y += self.y[0:n]

    def _refresh_properties(self):
        # self.properties.width = np.amax(np.sqrt([x ** 2 + y ** 2 for x, y in zip(self.x, self.y)]))
        self.properties.width = np.amax(np.sqrt(self.x ** 2 + self.y ** 2))
        self.properties.length = np.amax(self.z)
        self.properties.n_atoms = len(self.x)
        self.properties.density = len(self.x) / self.properties.length / (np.pi * self.properties.width ** 2)

        self.properties.width /= self.medium_atom.lbar
        self.properties.length /= self.medium_atom.lbar
        self.properties.density *= self.medium_atom.lbar ** 3
