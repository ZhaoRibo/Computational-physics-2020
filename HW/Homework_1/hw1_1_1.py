''' 正方形电阻网络 '''
from wheels import Matrix
from wheels import zero_mat
from wheels import up_tri_msolve
import time


def get_G(n):
    '''
    获取正方形电阻网络的导纳矩阵
    节点标号方式为从左下起，由左向右，依次取完每行
    return Matrix
    '''
    g = Matrix(zero_mat(n**2, n**2))
    for i in range(n**2):
        for j in range(n**2):
            if (abs(j - i) == 1) | (abs(j - i) == n):
                g[i][j] = -1
                g[i][i] += 1
    for i in range(1, n):
        g[i * n][i * n - 1] = 0
        g[i * n][i * n] -= 1
        g[i * n - 1][i * n] = 0
        g[i * n - 1][i * n - 1] -= 1
    return g


def solve_R(g, a, b, n):
    '''
    直接法求解标号为a,b节点间的等效电阻R
    return R
    '''
    assert a < b, "b>a不符合规范"

    # 删去g的第a行a列，使其成为非奇异矩阵
    g.pop_col(a)
    g.pop_row(a)

    m = g.row
    # 对称正定矩阵Cholesky分解改进算法，G=LDL^T,带状矩阵优化
    for i in range(m):
        if i < n:
            for j in range(i + 1):
                for k in range(j):
                    g[i][j] = g[i][j] - g[i][k] * g[j][k] / g[k][k]
        else:
            for j in range(i - n, i + 1):
                for k in range(j):
                    g[i][j] = g[i][j] - g[i][k] * g[j][k] / g[k][k]

    # 对电流I进行处理，I=LDI'
    # 求解 LDI'=I 得到 I'，下三角前代算法
    I = Matrix(zero_mat(m, 1))
    I[b - 1][0] = 1
    for i in range(m):
        assert g[i][i] != 0, "系数矩阵主元存在‘0’项"
        for j in range(i):
            I[i][0] = I[i][0] - g[i][j] / g[j][j] * I[j][0]
        I[i][0] = I[i][0]
    """
    # 利用(TD)^(-1)特性方法失败
    I = Matrix(zero_mat(m, 1))
    I[b-1][0] = 1
    for i in range(b, m):
        I[i][0] = -g[i][b - 1] / g[b - 1][b - 1]
    """

    g = g.T()
    # 上三角矩阵回代得解
   
    # 直接调用wheels通法
    x = up_tri_msolve(g, I)
    return x[b-1][0]
    """
    # 上三角回代带状优化, 直接覆盖在 I 上，存储优化
    for i in reversed(range(m)):
        assert g[i][i] != 0, "系数矩阵主元存在‘0’项"
        if i >= m - 1 - n:
            for j in range(m - 1, i, -1):
                I[i][0] = I[i][0] - g[i][j] * I[j][0]
        else:
            for j in range(i + n, i, -1):
                I[i][0] = I[i][0] - g[i][j] * I[j][0]
        I[i][0] = I[i][0] / g[i][i]

    return I[b - 1][0]
    """


start = time.time()
g = get_G(2)
print('边长: 1; ac')
print('R= ', solve_R(g, 0, 3, 2))
end = time.time()
print('datatime: ', end - start, '\n')

start = time.time()
g = get_G(2)
print('边长: 1; ab')
print('R= ', solve_R(g, 0, 1, 2))
end = time.time()
print('datatime: ', end - start, '\n')

start = time.time()
g = get_G(5)
print('边长: 4; ac')
print('R= ', solve_R(g, 0, 24, 5))
end = time.time()
print('datatime: ', end - start, '\n')

start = time.time()
g = get_G(5)
print('边长: 4; ab')
print('R= ', solve_R(g, 0, 4, 5))
end = time.time()
print('datatime: ', end - start, '\n')

start = time.time()
g = get_G(17)
print('边长: 16; ac')
print('R= ', solve_R(g, 0, 288, 17))
end = time.time()
print('datatime: ', end - start, '\n')

start = time.time()
g = get_G(17)
print('边长: 16; ab')
print('R= ', solve_R(g, 0, 16, 17))
end = time.time()
print('datatime: ', end - start, '\n')

start = time.time()
g = get_G(65)
print('边长: 64; ac')
print('R= ', solve_R(g, 0, 4224, 65))
end = time.time()
print('datatime: ', end - start, '\n')

start = time.time()
g = get_G(65)
print('边长: 64; ab')
print('R= ', solve_R(g, 0, 64, 65))
end = time.time()
print('datatime: ', end - start, '\n')
