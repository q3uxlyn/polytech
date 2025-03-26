import numpy as np
import matplotlib.pyplot as plt

# Заданная матрица 4x4
A = np.array([
    [2, 0, 0, -1],
    [0, 3, 0, -2],
    [1, -4, 6, 0],
    [0, -2, 2, 5]
])

# Размер матрицы
n = A.shape[0]

# Список цветов для кругов
colors = ['blue', 'green', 'red', 'purple']

# Построение кругов Гершгорина
fig, ax = plt.subplots()

for i in range(n):
    center = A[i, i]  # Центр круга (диагональный элемент)
    radius = np.sum(np.abs(A[i, :])) - np.abs(A[i, i])  # Радиус круга
    circle = plt.Circle((center.real, center.imag), radius, color=colors[i], fill=False)
    ax.add_patch(circle)
    plt.scatter(center.real, center.imag, color=colors[i])  # Центр круга

# Настройка графика
ax.set_aspect('equal', adjustable='datalim')
plt.xlim([-15, 15])  # Ограничения по оси X
plt.ylim([-10, 10])  # Ограничения по оси Y
plt.xlabel('Re')
plt.ylabel('Im')
plt.grid(True)
plt.show()
