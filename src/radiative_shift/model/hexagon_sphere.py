import numpy as np
from .hexagon import HexagonModel
from .param import SPHERE_RADIUS
from .param import RADIUS


class HexagonSphere(HexagonModel):

    def __init__(self, layers, copies, radius = SPHERE_RADIUS):
        # Просто вызывает конструктор цилиндра
        ncylinder = (copies + 1) * (1 + 3 * (layers + 2) * (layers + 1))
        rcylinder = (layers + 1) * RADIUS
        super().__init__(ncylinder, rcylinder)
        excess = []
        zm = np.average(self.z)
        for i, (xi, yi, zi) in enumerate(zip(self.x, self.y, self.z)):
            if xi ** 2 + yi ** 2 + (zi - zm) ** 2 > radius ** 2:
                excess.append(i)

        self.x = np.delete(self.x, excess)
        self.y = np.delete(self.y, excess)
        self.z = np.delete(self.z, excess)
