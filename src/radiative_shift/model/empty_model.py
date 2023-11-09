import numpy as np
from .general import GeneralModel
from .param import LBAR


class EmptyModel(GeneralModel):

    def __init__(self):
        super().__init__()
        self.measure_properties()

    def measure_properties(self):
        # TODO: write normal method
        # self.properties.radius = np.amax(np.sqrt([x ** 2 + y ** 2 for x, y in zip(self.x, self.y)])) / LBAR
        # self.properties.length = np.amax(self.z) / LBAR
        self.properties.radius = 1
        self.properties.length = 1
        self.properties.noa = len(self.x)
        # density = n0 * LBAR ** 3
        self.properties.density = len(self.x) / self.properties.length / (self.properties.radius ** 2)
