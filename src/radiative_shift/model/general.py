import numpy as np
from .param import LFACTOR
from abc import ABC
import logging


class Properties:
    radius: float
    length: float
    density: float
    noa: float


class GeneralModel(ABC):
    x: [].__class__
    y: [].__class__
    z: [].__class__
    properties: Properties

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


    def measureProperties(self):

        self.properties = Properties()
        self.properties.radius = np.amax(np.sqrt([x**2 + y**2 for x, y in zip(self.x, self.y)])) / 2 / np.pi
        self.properties.length = np.amax(self.z)
        self.properties.noa = len(self.x)
        self.properties.density = len(self.x) / self.properties.length / self.properties.radius**2 / np.pi

    def writeLog(self):

        logging.info("=========================================================")
        logging.info("Сylinder parameters")
        logging.info("Length = {:2.2}".format(self.properties.length) + " λ / 2 π")
        logging.info("Radius = {:2.2}".format(self.properties.radius) + " λ / 2 π")
        logging.info("Number of atoms {:}".format(self.properties.noa))
        logging.info("Density {:2.2}".format(self.properties.density) + " nλbar^3")
        logging.info("=========================================================")
