import math


# При w < 1 нижняя релаксация, при w > 1 верхняя релаксация, при w = 1 метод Зейделя
def relaxation(matrix_a, vector_b, w, eps):
    n = len(vector_b)
    vector_x = [0.0] * n
    vector_xn = [0.0] * n

    k = 0
    norm = 0
    while not norm > eps:
        k += 1
        norm = 0

        for i in range(n):
            vector_x[i] = vector_b[i]
            for j in range(n):
                if i != j:
                    vector_x[i] = vector_x[i] - matrix_a[i][j] * vector_x[j]
            vector_x[i] /= matrix_a[i][i]

            vector_x[i] = w*vector_x[i]+(1-w)*vector_xn[i]

            if math.fabs(vector_x[i] - vector_xn[i]) > norm:
                norm = math.fabs(vector_x[i]-vector_xn[i])
            vector_xn[i] = vector_x[i]

    return vector_x
