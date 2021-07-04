import numpy as np
import random
from .param import RADIUS
from .param import LENGTH
from .general import GeneralModel


class DisorderedModel(GeneralModel):

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
            z = random.uniform(i * LENGTH * copies / n, (i + 1) * LENGTH * copies / n)
            self.x.append(x)
            self.y.append(y)
            self.z.append(z)
