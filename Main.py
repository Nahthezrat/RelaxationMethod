import math
import numpy as np
import RelaxationMethod
import MatrixOperations as mo

# Начальные данные
N = 13
matrix_a = [
    [10.0*N+1, 4.0, 2.0, 2.0],
    [4.0, 8.0, 0.0, 2.0],
    [2.0, 0.0, 9.0, -4.0],
    [2.0, 2.0, -4.0, 12.0]
]
vector_b = [
    2.0*N*math.sin(N),
    5.0*(math.sin(N)-math.cos(N)),
    7.0*(math.sin(N)+math.cos(N)),
    3.0*math.sin(N)]

# Выводим систему (эхо-вывод)
print("\nLinear System: ")
for row in range(len(vector_b)):
    print("(", end='')
    for col in range(len(matrix_a[row])):
        print("\t{1:10.4f}{0}".format(" ", matrix_a[row][col]), end='')
    print("\t) * (\tX{0}) = (\t{1:10.4f})".format(row+1, vector_b[row]))

# Решаем систему методом верхней релаксации
# При w < 1 нижняя релаксация, при w > 1 верхняя релаксация, при w = 1 метод Зейделя
vector_x = np.array(RelaxationMethod.relaxation(matrix_a, vector_b, 1.5, 0.00001))

# Выводим решение
print("\nSolution: ")
print("\n".join("X{0} =\t{1:10.4f}".format(i+1, x) for i, x in enumerate(vector_x)))

# Вычисляем невязку: разницу между точным решением и приближённым
# r = b - x * A
A = np.array(matrix_a)
B = np.array(vector_b)
vector_r = B - np.dot(vector_x, A)
print("\nResidual: ")
print("\n".join("R{0} =\t{1:10.4f}".format(i+1, x) for i, x in enumerate(vector_r)))

# Формат данных для вычисления детерминанта и построения обратной матрицы
data = []  # Создаем пустую матрицу
for i in range(len(matrix_a)):
    for j in range(len(matrix_a)):
        data = data+[matrix_a[i][j]]  # Переводим квадратную матрицу в плоскую
dim = len(matrix_a)
dim2 = dim * dim


# Определитель матрицы
print("\nDeterminant: ")
print(mo.determinant(dim, data))

# Обратная матрица
adj_matrix = mo.adjoint_matrix(dim, data)  # Находим значение строк матрицы
for n in range(dim2): # Создаем счетчик итерации
    cs = int((adj_matrix[n])/mo.determinant(dim, data)) # Находим значение ячейки матрицы

# Выдуваем обратную матрицу
adj_matrix_normal = [[0.0] * len(vector_b) for i in range(len(vector_b))]
for i in range(0, 15, 4):
    for j in range(0, 4):
        adj_matrix_normal[i//4][j] = adj_matrix[i+j]

# Выводим обратную матрицу
print("\nAdjoint Matrix: ")
for row in range(len(adj_matrix_normal)):
    print("(", end='')
    for col in range(len(adj_matrix_normal[row])):
        print("\t{1:10.4f}{0}".format(" ", adj_matrix_normal[row][col]), end='')
    print("\t)")

# Перемножаем с исходной для проверки
print("\nVerification (A * A(-1) = E: ")
print(np.dot(matrix_a, adj_matrix_normal))
