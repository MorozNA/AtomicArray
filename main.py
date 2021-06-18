import numpy as np
from src.radiative_shift import HexagonModel
from src.radiative_shift import markovianSigmaMatrixForV

# Теперь у нас есть доступ к нашим атомам здесь


model = HexagonModel(3, 5)

gammas, _ = np.linalg.eig(markovianSigmaMatrixForV(model))

# По идее здесь ошибка - не должно быть положительных значений, ошибка в markovian. Буду проверять.
print(np.sort(np.imag(gammas)))
