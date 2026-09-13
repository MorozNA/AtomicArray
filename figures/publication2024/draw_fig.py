from numpy import linspace, loadtxt, real, argmax, round, zeros, vdot
from draw import draw

x = linspace(0.01, 5.0, 200)
print(x[52])
y5 = loadtxt('data/fig5.txt', dtype=complex)
y6 = loadtxt('data/fig6.txt', dtype=complex)
y7 = loadtxt('data/fig7.txt', dtype=complex)
y8 = loadtxt('data/fig8.txt', dtype=complex)

eigv5 = loadtxt('data/fig5_eigv40.txt', dtype=complex)

prod = zeros((11, 11))
for i in range(11):
    print(f'eigenvalue {i}: ', round(eigv5[:, i], 2))
    for j in range(11):
        prod[i, j] = abs(vdot(eigv5[:, i], eigv5[:, j]))
print(round(prod, 2))

eigv6 = loadtxt('data/fig6_eigv.txt', dtype=complex)
eigv7 = loadtxt('data/fig7_eigv.txt', dtype=complex)
eigv8 = loadtxt('data/fig8_eigv.txt', dtype=complex)
ind5 = []
ind6 = []
ind7 = []
ind8 = []
for i in range(y5.shape[0]):
    ind5.append(argmax(eigv5[:, i]) - 5)
    ind6.append(argmax(eigv6[:, i]) - 5)
for i in range(y7.shape[0]):
    ind7.append(argmax(eigv7[:, i]) - 3)
    ind8.append(argmax(eigv8[:, i]) - 3)

print(ind5)
print(ind6)
print(ind7)
print(ind8)

monte_yCS = real(loadtxt('data/monte-carlo_cs.txt', dtype=complex))
monte_yRB = real(loadtxt('data/monte-carlo_rb.txt', dtype=complex))


x_monte = x
draw(x, y5, 12, 'Figure 5', ind5, 0, x_monte, monte_yCS)
draw(x, y6, 10, 'Figure 6', ind6, 0, x, monte_yCS)
draw(x, y7, 10, 'Figure 7', ind7, 0, x, monte_yRB)
draw(x, y8, 10, 'Figure 8', ind8, 2, x, monte_yRB)