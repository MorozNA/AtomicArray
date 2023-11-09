import numpy as np
from .general import GeneralModel
from .param import LBAR


# Конструктор класса строит заданное количество точек, равномерно распределенных внутри сферы радиуса radius.
class DisorderedSphere(GeneralModel):

    def __init__(self, density, radius):
        super().__init__()
        n = int(density * 4 / 3 * np.pi * radius ** 3)
        u = np.random.uniform(-1, 1, size=n)
        phi = np.random.uniform(0, 2 * np.pi, size=n)
        r = (np.random.uniform(0, 1, size=n) ** (1 / 3)) * radius
        self.x = r * np.cos(phi) * (1 - u ** 2) ** 0.5
        self.y = r * np.sin(phi) * (1 - u ** 2) ** 0.5
        self.z = r * u

        self.measure_properties()
        self.write_log()

    def measure_properties(self):
        self.properties.radius = np.amax(np.sqrt([x ** 2 + y ** 2 for x, y in zip(self.x, self.y)]))  # / LBAR
        self.properties.length = self.properties.radius
        self.properties.noa = len(self.x)
        self.properties.density = self.properties.noa / (4 / 3 * np.pi * self.properties.radius ** 3)
