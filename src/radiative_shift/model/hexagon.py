import numpy as np
from .general import GeneralModel
from .param import LBAR


class HexagonModel(GeneralModel):

    def __init__(self, length, radius, density):
        super().__init__()
        n = round(length * (density * np.pi * radius ** 2))
        self.unitradius = density ** (-1 / 3)
        self.unitlength = density ** (-1 / 3)
        layers = round(radius / self.unitradius)
        copies = round(length / self.unitlength)
        while True:
            if n >= copies * (1 + 3 * (layers + 1) * (layers + 2)):
                layers += 1
                self.unitradius = radius / layers
                if n >= (copies + 1) * (1 + 3 * (layers + 1) * layers):
                    copies += 1
                    self.unitlength = length / copies
                else:
                    break
            else:
                break
        self.unitradius = radius / layers
        self.unitlength = length / copies
        self._nl = 0
        # Lists are faster to append, after the addition of layers and copies they will be transformed to np.arrays
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

        self.measure_properties()
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
        self.z += n * [self.z[-1] + self.unitlength]
        self.x += self.x[0:n]
        self.y += self.y[0:n]

    def measure_properties(self):
        self.properties.radius = np.amax(np.sqrt([x ** 2 + y ** 2 for x, y in zip(self.x, self.y)])) / LBAR
        self.properties.length = np.amax(self.z) / LBAR
        self.properties.noa = len(self.x)
        self.properties.density = len(self.x) / self.properties.length / (np.pi * self.properties.radius ** 2)