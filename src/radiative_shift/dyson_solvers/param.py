import numpy as np

KV = 1
DDI = 1
HBAR = 1
GAMMA = 1


D01 = np.sqrt(KV**3 * HBAR*GAMMA/4)
D10 = D01
D00 = np.sqrt(KV**3 * HBAR*GAMMA/4)
D01M = np.sqrt(KV**3 * HBAR*GAMMA/4)
D1M0 = D01M