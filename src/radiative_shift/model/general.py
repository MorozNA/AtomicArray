import numpy as np


class GeneralModel:
    """
    TODO: юниттест
    """

    def __init__(self):
        self.x = []
        self.y = []
        self.z = []

    def getDistances(self):
        distanceXM = []  # Добавил матрицы в таком виде. Если нужно - поменяю на x +- iy
        distanceZ = []
        for i in range(len(self.x)):
            distanceXM.append([])
            distanceZ.append([])
            for j in range(len(self.x)):
                distanceXM[i].append((-self.x[i] + self.x[j] - 1j * self.y[i] + 1j * self.y[j]) / np.sqrt(2))
                distanceZ[i].append(self.z[i] - self.z[j])
        arrXM = np.asarray(distanceXM)
        arrZ = np.asarray(distanceZ)
        # Проверить
        arrR = np.sqrt(2 * np.square(np.abs(arrXM)) + np.square(arrZ)) \
               + np.identity(np.shape(arrZ)[0])  # Ignoring division by zero
        return arrXM / arrR, arrZ / arrR, arrR
