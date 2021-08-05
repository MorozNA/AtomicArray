import numpy as np
from .param import RADIUS
from .param import LENGTH
from .general import GeneralModel


class HexagonModel(GeneralModel):

    def __init__(self, n, radius, density):
        length = n / (density * np.pi * radius ** 2)
        self.unitradius = RADIUS
        self.unitlength = LENGTH
        layers = int(radius / self.unitradius)
        copies = int(length / self.unitlength)
        # Число атомов в одной копии: 1 + 3 * layers * (layers + 1)
        flag = 0
        while flag != 1:
            if n >= copies * (1 + 3 * (layers + 1) * (layers + 2)):
                layers += 1
                self.unitradius = radius / layers
                if n >= (copies + 1) * (1 + 3 * (layers + 1) * layers):
                    copies += 1
                    self.unitlength = length / copies
                else:
                    flag = 1
            else:
                flag = 1
        self._nl = 0  # n - число слоев. Меняется при добавлении слоя.
        self.x = [0]
        self.y = [0]
        self.z = [0]

        for _ in range(layers):
            self._addLayer()

        for _ in range(copies):
            self._addCopy()

        self.measureProperties()
        self.writeLog()

    def _addLayer(self):
        self._nl = self._nl + 1
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

    def _addCopy(self):
        n = int(1 + self._nl * (6 + 6 * self._nl) / 2)  # Тут считаю число атомов в одной плоскости
        self.z = self.z + n * [self.z[-1] + self.unitlength]
        self.x = self.x + self.x[0:n]
        self.y = self.y + self.y[0:n]