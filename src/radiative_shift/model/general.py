import numpy as np


class GeneralModel:

    def __init__(self):
        self.x = []
        self.y = []
        self.z = []

    def getDistances(self):
        distanceX = []  # Добавил матрицы в таком виде. Если нужно - поменяю на x +- iy
        distanceY = []
        distanceZ = []
        for i in range(len(self.x)):
            distanceX.append([])
            distanceY.append([])
            distanceZ.append([])
            for j in range(len(self.x)):
                distanceX[i].append(self.x[i] - self.x[j])
                distanceY[i].append(self.y[i] - self.y[j])
                distanceZ[i].append(self.z[i] - self.z[j])
        arrX = np.asarray(distanceX)
        arrY = np.asarray(distanceY)
        arrZ = np.asarray(distanceZ)
        arrR = np.sqrt(np.square(arrX) + np.square(arrY) + np.square(arrZ))
        return arrX, arrY, arrZ, arrR
