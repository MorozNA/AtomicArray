import numpy as np
from src.radiative_shift import HexagonModel, HexagonSphere
from src.radiative_shift import markovianSigmaMatrixForV

# Теперь у нас есть доступ к нашим атомам здесь

model = HexagonModel(5, 5)

gammas, _ = np.linalg.eig(markovianSigmaMatrixForV(model))

print(np.sort(np.imag(gammas)))
