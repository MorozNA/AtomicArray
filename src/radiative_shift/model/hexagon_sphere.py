import numpy as np
from .hexagon import HexagonModel
from .param import SPHERE_RADIUS


class HexagonSphere(HexagonModel):

    def __init__(self, layers, copies):
        # Просто вызывает конструктор цилиндра
        super().__init__(layers, copies)
        excess = []
        zm = np.average(self.z)
        for i, (xi, yi, zi) in enumerate(zip(self.x, self.y, self.z)):
            if xi ** 2 + yi ** 2 + (zi - zm) ** 2 > SPHERE_RADIUS ** 2:
                excess.append(i)

        self.x = np.delete(self.x, excess)
        self.y = np.delete(self.y, excess)
        self.z = np.delete(self.z, excess)
