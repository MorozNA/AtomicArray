import numpy as np
from .general import GeneralModel
from .param import LBAR


class DisorderedModel(GeneralModel):

    def __init__(self, length, radius, density):
        super().__init__()
        n = int(length * (density * np.pi * radius ** 2))

        # Generate random points using numpy's random functions
        phi = np.random.uniform(0, 2 * np.pi, size=n)
        u = np.random.uniform(0, radius ** 2, size=n)
        x = np.sqrt(u) * np.cos(phi)
        y = np.sqrt(u) * np.sin(phi)
        z = np.random.uniform(0, length, size=n)

        # TODO: x, y and z are already np.arrays
        # Use numpy arrays instead of Python lists
        self.x = np.array(x)
        self.y = np.array(y)
        self.z = np.array(z)

        self.measure_properties()

    def measure_properties(self):
        self.properties.radius = np.amax(np.sqrt([x ** 2 + y ** 2 for x, y in zip(self.x, self.y)])) / LBAR
        self.properties.length = np.amax(self.z) / LBAR
        self.properties.noa = len(self.x)
        self.properties.density = len(self.x) / self.properties.length / (np.pi * self.properties.radius ** 2)
