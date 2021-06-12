import numpy as np
import matplotlib.pyplot as plt

radius = 1
length = 1


class AtomicArray:
    def __init__(self, n):
        self.n = n + 1  # n - число точек в одной плоскости. Оно меняется при добавлении слоя.
        self.angle = 2 * np.pi / n  # angle - определяет какой именно многоуголник мы повторяем. Не меняется.
        self.x = [0]
        self.y = [0]
        for v in range(n):
            x = radius * np.sin(v * self.angle)  # В цикле создаём сам многоугольник вокруг начала координат.
            y = radius * np.cos(v * self.angle)
            self.x.append(x)
            self.y.append(y)
        self.z = [0] * (n + 1)

    def getX(self):
        return self.x

    def getY(self):
        return self.y

    def getZ(self):
        return self.z

    def addLayer(self):  # Пытался улучшить этот метод. Пока не знаю как. Не нравится, что проход по всем точкам идёт.
        flag = 0
        for i in range(1, len(self.x)):
            for j in range(len(self.x)):
                x = self.x[i] + radius * np.sin(j * self.angle)
                y = self.y[i] + radius * np.cos(j * self.angle)
                for v in range(len(self.x)):
                    if self.x[v] + 0.01 > x > self.x[v] - 0.01 and self.y[v] + 0.01 > y > self.y[v] - 0.01 and self.z[
                        v] == self.z[i]:
                        flag = 0
                        break
                    else:
                        flag = 1
                if flag == 1:
                    self.x.append(x)
                    self.y.append(y)
                    self.z.append(self.z[i])
                    self.n = self.n + 1
                else:
                    pass

    def addCopy(self):
        self.z = self.z + self.n * [self.z[-1] + 1]
        self.x = self.x + self.x[0:self.n]
        self.y = self.y + self.y[0:self.n]


# Тест программы. Можно сначала добавлять копии, но тогда addLayer будет проходить по всем копиям, что занимает время.
if __name__ == '__main__':
    test = AtomicArray(6)
    layers = 8
    copies = 18

    for i in range(layers):
        test.addLayer()

    for j in range(copies):
        test.addCopy()

    x = test.getX()
    y = test.getY()
    z = test.getZ()
    print(len(x), len(y), len(z))  # Проверял не добавил ли случайно лишних точек. Заменил '== x' в addLayer на '> x >'

    plt.plot(x, y, 'o')
    plt.show()

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.scatter3D(x, y, z)
    plt.show()
