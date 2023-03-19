import numpy as np
from .general import GeneralModel


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

        # Use numpy arrays instead of Python lists
        self.x = np.array(x)
        self.y = np.array(y)
        self.z = np.array(z)

        self.measure_properties()
        self.write_log()