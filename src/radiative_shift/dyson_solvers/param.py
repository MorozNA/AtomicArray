import numpy as np

N_refr = 1.45

DDI = 1
HBAR = 6.62607015e-27 / 2 / np.pi
C = 2.998e10

# Rubidium 87
LBAR = 780e-7 / 2 / np.pi
OM = C / LBAR
KV = OM / C
GAMMA = 38.11e6 / 3

# Cesium 137
# LBAR = 852e-7 / 2 / np.pi
# OM = C / LBAR
# KV = OM / C
# GAMMA = 32.815e6

# Rubidium 87
# GAMMA = 1
# OM = (C / (780e-7 / 2 / np.pi)) / (38.11e6) * GAMMA
# LBAR = C / OM
# KV = OM / C