import numpy as np
from .hexagon import HexagonModel


class HexagonSphere(HexagonModel):

    def __init__(self, density, radius):
        super().__init__(2 * radius, radius, density)
        excess = []
        zm = np.average(self.z)
        for i, (xi, yi, zi) in enumerate(zip(self.x, self.y, self.z)):
            if xi ** 2 + yi ** 2 + (zi - zm) ** 2 > radius ** 2:
                excess.append(i)

        self.x = np.delete(self.x, excess)
        self.y = np.delete(self.y, excess)
        self.z = np.delete(self.z, excess)
