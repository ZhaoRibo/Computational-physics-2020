''' 正方形电阻网络 '''
from wheels import Matrix
from wheels import zero_mat
import time
import copy


def get_G(n):
    '''
    获取六边形电阻网络的导纳矩阵
    进行星->三角变换后变为等边长三角形网络
    取边上节点数为n，节点标号方式为从一个顶点开始，算作第一层，一次向上取第二层，第三层……
    return Matrix
    '''
    g = Matrix(zero_mat(int((n ** 2 + n) / 2), int((n ** 2 + n) / 2)))
    g[0][0] = 1
    g[0][1] = -1/2
    g[0][2] = -1/2
    i = 1
    for k in range(1, n-1):
        g[i][i - k] = -1/2
        g[i][i + 1] = -1/3
        g[i][i + k + 1] = -1/2
        g[i][i + k + 2] = -1/3
        g[i][i] = 5/3
        for j in range(1, k):
            s = i + j
            g[s][s - k - 1] = -1/3
            g[s][s - k] = -1/3
            g[s][s - 1] = -1/3
            g[s][s + 1] = -1/3
            g[s][s + k + 1] = -1/3
            g[s][s + k + 2] = -1/3
            g[s][s] = 2
        s = i + k
        g[s][s - k - 1] = -1/2
        g[s][s - 1] = -1/3
        g[s][s + k + 1] = -1/3
        g[s][s + k + 2] = -1/2
        g[s][s] = 5/3
        i = s + 1
    g[i][i - n + 1] = -1/2
    g[i][i + 1] = -1/2
    g[i][i] = 1
    for j in range(1, n - 1):
        s = i + j
        g[s][s - n] = -1/3
        g[s][s - n + 1] = -1/3
        g[s][s - 1] = -1/2
        g[s][s + 1] = -1/2
        g[s][s] = 5/3
    s = i + n - 1
    g[s][s - n] = -1/2
    g[s][s - 1] = -1/2
    g[s][s] = 1
    # 带状矩阵，半宽为n
    return g


def col_abs(x):
    '''
    x: Martix; x.col=1
    列向量x的欧几里得范数
    '''
    return (x.T() * x)[0][0]


def solve_R(g, a, b, n, M, e):
    '''
    迭代法求解标号为a,b节点间的等效电阻R
    return R
    '''
    assert a < b, "b>a不符合规范"

    # 删去g的第a行a列，使其成为非奇异矩阵
    g.pop_col(a)
    g.pop_row(a)

    m = g.row

    # 对电流I进行处理
    I = Matrix(zero_mat(m, 1))
    I[b - 1][0] = 1
    """
    # Jacobi迭代法收敛速度太慢
    # Jacobi迭代法求解 gx=I
    x_0 = Matrix(zero_mat(m, 1))
    for i in range(m):
        x_0[i][0] = 1           # 内置初始解向量
    x = Matrix(zero_mat(m, 1))  # 另一组解向量初始化
    M = 1000     # 内置最大迭代次数100
    e = 1e-10  # 内置判停标准1e-10
    for k in range(M):
        for i in range(m):
            x[i][0] = I[i][0]
            for j in range(m):
                x[i][0] = x[i][0] - g[i][j] * x_0[j][0]
            x[i][0] = x[i][0] / g[i][i] + x_0[i][0]
        if col_abs(x - x_0) < e:
            print("迭代次数: ", k+1)
            break
        elif k == M - 1:
            print("Not converged at given M=", M, " and e=", e)
            break
        x_0 = copy.deepcopy(x)
    """
    # Gauss-Seidel迭代法 + 带状对称矩阵优化
    x = Matrix(zero_mat(m, 1))
    for i in range(m):
        x[i][0] = 1  # 内置初始解向量
    x_0 = Matrix(zero_mat(m, 1))  # 另一组解向量初始化
    for k in range(M):
        x_0 = copy.deepcopy(x)
        if m > 2 * n + 1:
            for i in range(m):
                x[i][0] = I[i][0]
                if i < n:
                    for j in range(i):
                        x[i][0] = x[i][0] - g[i][j] * x[j][0]
                    for j in range(i + 1, i + n + 1):
                        x[i][0] = x[i][0] - g[i][j] * x[j][0]
                    x[i][0] = x[i][0] / g[i][i]
                elif i > m - n - 1:
                    for j in range(i - n, i):
                        x[i][0] = x[i][0] - g[i][j] * x[j][0]
                    for j in range(i + 1, m):
                        x[i][0] = x[i][0] - g[i][j] * x[j][0]
                    x[i][0] = x[i][0] / g[i][i]
                else:
                    for j in range(i - n, i):
                        x[i][0] = x[i][0] - g[i][j] * x[j][0]
                    for j in range(i + 1, i + n + 1):
                        x[i][0] = x[i][0] - g[i][j] * x[j][0]
                    x[i][0] = x[i][0] / g[i][i]
                # 逐次超松弛迭代法
                x[i][0] = -0.9 * x_0[i][0] + 1.9 * x[i][0]
        else:
            for i in range(m):
                x[i][0] = I[i][0]
                for j in range(i):
                    x[i][0] = x[i][0] - g[i][j] * x[j][0]
                for j in range(i + 1, m):
                    x[i][0] = x[i][0] - g[i][j] * x[j][0]
                x[i][0] = x[i][0] / g[i][i]
                # 逐次超松弛迭代法
                x[i][0] = -0.9 * x_0[i][0] + 1.9 * x[i][0]
        if col_abs(x - x_0) < e:
            print("迭代次数: ", k + 1)
            break
        elif k == M - 1:
            print("Not converged at given M=", M, " and e=", e)
            break

    return x[b - 1][0]


start = time.time()
g = get_G(5)
print('边长: 4; ac')
print('R= ', solve_R(g, 0, 9, 5, 1000, 1e-10))
end = time.time()
print('datatime: ', end - start, '\n')


start = time.time()
g = get_G(17)
print('边长: 16; ac')
print('R= ', solve_R(g, 0, 135, 17, 1000, 1e-10))
end = time.time()
print('datatime: ', end - start, '\n')


start = time.time()
g = get_G(65)
print('边长: 64; ac')
print('R= ', solve_R(g, 0, 2079, 65, 2000, 1e-3))
end = time.time()
print('datatime: ', end - start, '\n')
