import numpy as np
from .general import GeneralModel


class HexagonModel(GeneralModel):

    def __init__(self, length, radius, density):
        super().__init__()
        n = int(length * (density * np.pi * radius ** 2))
        self.unitradius = (radius / int(radius / (density ** (1/3)))).real
        self.unitlength = (length / int(length / (density ** (1/3)))).real
        layers = int(radius / self.unitradius)
        copies = int(length / self.unitlength)
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
        self._nl = 0
        self.x = [0]
        self.y = [0]
        self.z = [0]
        for _ in range(layers):
            self._add_layer()
        for _ in range(copies):
            self._add_copy()
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
        n = int(1 + self._nl * (6 + 6 * self._nl) / 2)
        self.z += n * [self.z[-1] + self.unitlength]
        self.x += self.x[0:n]
        self.y += self.y[0:n]
