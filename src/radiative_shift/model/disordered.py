import numpy as np
from .param import RADIUS
from .param import LENGTH
import random


class DisorderedModel:

    def __init__(self, layers, copies):
        # число точек в одной плоскости (1 + 3 * (layers + 2) * (layers + 1))
        # число плоскостей (1 + copies)
        n = (copies + 1) * (1 + 3 * (layers + 2) * (layers + 1))
        self.x = []
        self.y = []
        self.z = []
        for i in range(n):
            phi = random.uniform(0, 2 * np.pi)
            u = random.uniform(0, ((layers + 1) * RADIUS) ** 2)
            x = np.sqrt(u) * np.cos(phi)
            y = np.sqrt(u) * np.sin(phi)
            z = random.uniform(0, LENGTH * copies)
            self.x.append(x)
            self.y.append(y)
            self.z.append(z)

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
