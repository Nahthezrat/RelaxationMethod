import math


# Вычисление определителя
def determinant(dim, data):
    if dim == 0:
        det = 0
    if dim == 1:
        det = data[0]
    if dim == 2:
        det = data[0] * data[3] - data[2] * data[1]
    if dim == 3:
        det = data[0] * data[4] * data[8] + data[1] * data[5] * data[6] + data[2] * data[3] * data[7] - data[6] * data[4] * data[2] - data[7] * data[5] * data[0] - data[8] * data[3] * data[1]
    if dim > 3:
        det = 0
        i = 0
        for j in range(dim):
            det += data[i * dim + j] * cofactor(dim, data, i, j)
    return det


# Вычисление минора
def minor(dim, data, i, j): 
    new_matrix = []
    for  k in range(dim*dim):
        if math.floor(k / dim) != i and k % dim != j:
            new_matrix.append(data[k])
    return determinant(dim - 1, new_matrix)


# Вычисление кофактора
def cofactor(dim, data, i, j):
    return math.pow(-1, i + j) * minor(dim, data, i, j)


# Вычисление кофакторной матрицу
def cofactor_matrix(dim, data):
    new_matrix = []
    for  i in range(dim):
        for  j in range(dim):
            new_matrix.append(cofactor(dim, data, i, j))
    return new_matrix


# Транспонирование матрицы
def transpose_matrix(dim, data):
    new_matrix = []
    for i in range(dim):
        for j in range(dim):
            new_matrix.append(data[j * dim + i])
    return new_matrix


# Вычисление обратной матрицы
def adjoint_matrix(dim, data):
    return transpose_matrix(dim, cofactor_matrix(dim, data))
