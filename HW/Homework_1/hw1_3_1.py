''' 六边形电阻网络 '''
from wheels import Matrix
from wheels import zero_mat
from wheels import up_tri_msolve
import time


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
g = get_G(5)
print('边长: 4; ac')
print('R= ', solve_R(g, 0, 9, 5))
end = time.time()
print('datatime: ', end - start, '\n')


start = time.time()
g = get_G(17)
print('边长: 16; ac')
print('R= ', solve_R(g, 0, 135, 17))
end = time.time()
print('datatime: ', end - start, '\n')


start = time.time()
g = get_G(65)
print('边长: 64; ac')
print('R= ', solve_R(g, 0, 2079, 65))
end = time.time()
print('datatime: ', end - start, '\n')