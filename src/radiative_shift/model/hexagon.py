import numpy as np
from .param import RADIUS
from .param import LENGTH


class HexagonModel:

    def __init__(self, layers, copies):
        self._nl = 1  # n - число слоев. Меняется при добавлении слоя.
        self.x = [0]
        self.y = [0]
        for v in range(6):
            x = RADIUS * np.sin(v * 2 * np.pi / 6)  # Везде убрал self._angle, не думаю что понадобится эта переменная
            y = RADIUS * np.cos(v * 2 * np.pi / 6)
            self.x.append(x)
            self.y.append(y)
        self.z = [0] * (6 + 1)

        for _ in range(layers):
            self._addLayer()

        for _ in range(copies):
            self._addCopy()

        self.distanceX = [[]]  # Добавил матрицы в таком виде. Если нужно - поменяю на x +- iy
        self.distanceY = [[]]
        self.distanceZ = [[]]
        self.distanceR = [[]]
        for i in range(len(self.x)):
            for j in range(len(self.x)):
                self.distanceX[i].append(self.x[i] - self.x[j])
                self.distanceY[i].append(self.y[i] - self.y[j])
                self.distanceZ[i].append(self.z[i] - self.z[j])
                self.distanceR[i].append(
                    np.sqrt(self.distanceX[i][j] ** 2 + self.distanceY[i][j] ** 2 + self.distanceZ[i][j] ** 2))
            self.distanceX.append([])
            self.distanceY.append([])
            self.distanceZ.append([])
            self.distanceR.append([])

    def _addLayer(self):  # Получилось улучшить метод, теперь нет прохода по точкам. Новые точки просто добавляются.
        self._nl = self._nl + 1
        layerX = []
        layerY = []
        layerZ = []
        for i in range(6):
            x = self._nl * RADIUS * np.sin(i * 2 * np.pi / 6)
            y = self._nl * RADIUS * np.cos(i * 2 * np.pi / 6)
            z = 0
            layerX.append(x)
            layerY.append(y)
            layerZ.append(z)
            for j in range(1, self._nl):
                layerX.append(x + j * RADIUS * np.sin((i + 2) * 2 * np.pi / 6))
                layerY.append(y + j * RADIUS * np.cos((i + 2) * 2 * np.pi / 6))
                layerZ.append(z)
        self.x += layerX
        self.y += layerY
        self.z += layerZ

    def _addCopy(self):
        n = int(1 + self._nl * (6 + 6 * self._nl) / 2)  # Тут считаю число атомов в одной плоскости
        self.z = self.z + n * [self.z[-1] + LENGTH]
        self.x = self.x + self.x[0:n]
        self.y = self.y + self.y[0:n]
