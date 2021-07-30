import numpy as np
from .param import LFACTOR
from abc import ABC


class GeneralModel(ABC):
    x: [].__class__
    y: [].__class__
    z: [].__class__

    def getDistances(self):
        distanceXM = []
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

    def removeDuplicates(self, r=LFACTOR):
        excess = []
        for i in range(len(self.x)):
            for j in range(len(self.x)):
                if i != j and (self.x[i] - self.x[j]) ** 2 \
                        + (self.y[i] - self.y[j]) ** 2 + (self.z[i] - self.z[j]) ** 2 < r ** 2:
                    excess.append(i)
                    continue
        self.x = np.delete(self.x, excess)
        self.y = np.delete(self.y, excess)
        self.z = np.delete(self.z, excess)
