import numpy as np
from .general import GeneralModel


class DisorderedComb(GeneralModel):

    # TODO: parameters should be changed to -> number of etches, period, density
    def __init__(self, length, period, density):
        super().__init__()
        a = period
        height = 0.5 * a
        width = 0.5 * a
        length_etched = 1.5 * a
        n = round(length * height * width * density)

        # Generate random points using numpy's random functions
        x = np.random.uniform(0, width, size=n).tolist()
        y = np.random.uniform(-height / 2, height / 2, size=n).tolist()
        z = np.random.uniform(0, length, size=n).tolist()

        # Create etched regions on one side of the cuboid
        num_etched = round(length / period)
        n_etched = round(length_etched * height * width * density)
        shift = (period - width) / 2
        # TODO: add x and y point;
        # TODO: each etch should have another number of points (it's volume * density)
        # TODO: change +=***.tolist() to np.append()
        for i in range(num_etched):
            z += np.random.uniform(shift + i * period, shift + i * period + width, size=n_etched).tolist()
            x += np.random.uniform(width, length_etched, size=n_etched).tolist()
            y += np.random.uniform(-height / 2, height / 2, size=n_etched).tolist()

        # Use numpy arrays instead of Python lists
        self.x = -np.array(x)
        self.y = np.array(y)
        self.z = np.array(z)

        self.measure_properties()
        self.write_log()

    def measure_properties(self):
        # TODO: maybe add another property: height
        self.properties.length = np.amax(self.z)  # / LBAR
        self.properties.radius = np.amax(abs(self.x))
        self.properties.noa = len(self.x)
        self.properties.density = len(self.x) / self.properties.length / self.properties.radius ** 2 / np.pi