import numpy as np
import random
from .param import RADIUS
from .param import LENGTH
from .general import GeneralModel


class DisorderedModel(GeneralModel):

    def __init__(self, n, radius):
        unitradius = radius / int(radius / RADIUS)
        unitlength = LENGTH
        layers = int(radius / unitradius - 1)
        copies = int(n / (1 + 3 * (layers + 2) * (layers + 1)))
        length = unitlength * copies
        # Если вводить density, то length = n / (density * np.pi * r ** 2)
        self.x = []
        self.y = []
        self.z = []
        for i in range(n):
            phi = random.uniform(0, 2 * np.pi)
            u = random.uniform(0, radius ** 2)
            x = np.sqrt(u) * np.cos(phi)
            y = np.sqrt(u) * np.sin(phi)
            z = random.uniform(0, length)
            self.x.append(x)
            self.y.append(y)
            self.z.append(z)
