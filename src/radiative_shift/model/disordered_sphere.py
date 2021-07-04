import numpy as np
import random
from .param import SPHERE_RADIUS
from .general import GeneralModel


class DisorderedSphere(GeneralModel):

    def __init__(self, n, radius=SPHERE_RADIUS):
        # число точек в одной плоскости (1 + 3 * (layers + 2) * (layers + 1))
        # число плоскостей (1 + copies)
        self.x = []
        self.y = []
        self.z = []
        for i in range(n):
            u = 2 * random.uniform(0, 1) - 1
            phi = 2 * np.pi * random.uniform(0, 1)
            r = random.uniform(0, radius) ** (1 / 3.)
            x = r * np.cos(phi) * (1 - u ** 2) ** 0.5
            y = r * np.sin(phi) * (1 - u ** 2) ** 0.5
            z = r * u
            self.x.append(x)
            self.y.append(y)
            self.z.append(z)
