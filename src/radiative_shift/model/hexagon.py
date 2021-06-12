import numpy as np
from .param import RADIUS
from .param import LENGTH


class HexagonModel:

    """
    TODO: нужны матрицы взаимных расстояний, взаимных разностей z, x + iy, x - iy (4 матрицы NxN)
    """
    def __init__(self, layers, copies):
        self._n = 7  # n - число точек в одной плоскости. Оно меняется при добавлении слоя.
        self._angle = 2 * np.pi / 6  # angle - определяет какой именно многоуголник мы повторяем. Не меняется.
        self.x = [0]
        self.y = [0]
        for v in range(6):
            x = RADIUS * np.sin(v * self._angle)  # В цикле создаём сам многоугольник вокруг начала координат.
            y = RADIUS * np.cos(v * self._angle)
            self.x.append(x)
            self.y.append(y)
        self.z = [0] * (6 + 1)

        for _ in range(layers):
            self._addLayer()

        for _ in range(copies):
            self._addCopy()

    # геттеры и сеттеры не нужны тут
    # если хочется как в enterprise, то нужно пометить поля чертой self.x -> self._x,
    # что будет значить, что это поле относится к приватным

    def _addLayer(self):  # Пытался улучшить этот метод. Пока не знаю как. Не нравится, что проход по всем точкам идёт.
        layerX = []
        layerY = []
        layerZ = []
        for i, _ in enumerate(self.x):
            # Строим шестиугольник в каждой вершине, можно улучшить в 3 раза, но не нужно
            for j in range(6):
                x = self.x[i] + RADIUS * np.sin(j * self._angle)
                y = self.y[i] + RADIUS * np.cos(j * self._angle)
                z = 0

                if (round(x, 2), round(y, 2)) \
                        in [(round(x, 2), round(y, 2)) for x, y
                            in zip(self.x + layerX, self.y + layerY)]:
                    continue
                else:
                    layerX.append(x)
                    layerY.append(y)
                    layerZ.append(z)
                    self._n = self._n + 1

        self.x += layerX
        self.y += layerY
        self.z += layerZ

    def _addCopy(self):
        self.z = self.z + self._n * [self.z[-1] + LENGTH]
        self.x = self.x + self.x[0:self._n]
        self.y = self.y + self.y[0:self._n]