import numpy as np
from .hexagon import HexagonModel
from .param import LBAR


class HexagonSphere(HexagonModel):

    def __init__(self, radius, density):
        super().__init__(2 * radius, radius, density)
        excess = []
        zm = np.average(self.z)
        for i, (xi, yi, zi) in enumerate(zip(self.x, self.y, self.z)):
            if xi ** 2 + yi ** 2 + (zi - zm) ** 2 > radius ** 2:
                excess.append(i)

        self.x = np.delete(self.x, excess)
        self.y = np.delete(self.y, excess)
        self.z = np.delete(self.z, excess)

        self.measure_properties()
        self.write_log()

    def measure_properties(self):
        self.properties.radius = np.amax(np.sqrt([x ** 2 + y ** 2 for x, y in zip(self.x, self.y)])) / LBAR
        self.properties.length = self.properties.radius
        self.properties.noa = len(self.x)
        self.properties.density = len(self.x) / (4 / 3 * np.pi * self.properties.radius ** 3)
