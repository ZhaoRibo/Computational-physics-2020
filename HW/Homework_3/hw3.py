import matplotlib.pyplot as plt
import copy
from wheels import Matrix, zero_mat
from math import pi, sin


def coefficientmatrix(n):
    """
    系数矩阵生成, 已包含边界条件‘u=0’
    input： 剖分数  n
    output：系数矩阵A
    """
    n0 = (n - 1) * (n - 1)              # 系数矩阵的维数
    A = Matrix(zero_mat(n0, n0))
    for i in range(0, n0):              # 设置中央对角线
        A[i][i] = 4
    for j in range(0, n - 1):           # 设置里面的两条斜对角线
        for k in range(0, n - 2):
            kk = j * (n - 1) + k
            A[kk][kk + 1] = -1
            A[kk + 1][kk] = -1
    for k in range(0, n0 - n + 1):      # 设置外面的两条斜对角线
        A[k][k + n - 1] = -1
        A[k + n - 1][k] = -1
    return A


def fxy(x, y):
    """泊松方程右侧非齐次项"""
    f2 = 2. * pi * pi * sin(pi * x) * sin(pi * y)
    return f2


def fij(n, h, x0, y0):
    """
    线性方程组右侧列向量
    input：剖分数n，步长h，x0，y0
    output：Matrix fij
    """
    n0 = (n - 1) * (n - 1)
    f = Matrix(zero_mat(n0, 1))
    for i in range(1, n):
        for j in range(1, n):
            k = (i - 1) * (n - 1) + j - 1
            f[k][0] = fxy(x0 + i * h, y0 + j * h) * h * h
    return f


def cg(A, b, x):
    """
    共轭梯度法求解线性方程组
    input：系数矩阵A，列向量b，初始值x
    output：解向量x
    """
    assert isinstance(A, Matrix) | isinstance(b, Matrix) | isinstance(x, Matrix), "输入数据类型错误"
    r = b - A * x           # r=b-Ax，r是梯度方向
    p = copy.deepcopy(r)    # 搜索方向
    i = 0                   # 记录迭代步数
    while ((r.T() * r)**0.5 > 1.e-10 and i < 100):
        a = (r.T() * p) / (r.T() * A * p)       # 搜索步长
        x = x + a * p                           # 更新解向量
        r = r - a * A * p                       # 更新残差向量/梯度向量
        b = -(r.T() * A * p) / (p.T() * A * p)
        p = r + b * p                           # 更新搜索方向
        i = i + 1
    return x


def e_l2(uij, n, h, x0, y0):
    """
    L2范数下的误差估计(高斯积分)
    input: 解矩阵 uij， 切分数 n， 步长 h， 初值 x0,y0
    output: 误差值 el2
    """
    el2 = 0
    for i in range(n):
        for j in range(n):
            el2 += singlecell(uij[i][j], uij[i + 1][j], uij[i + 1][j + 1], uij[i][j + 1], x0 + i * h, y0 + j * h)
    el2 = el2 ** 0.5 * h
    return el2


def singlecell(u1, u2, u3, u4, xi, yj):
    """
    一个单元格对 (uij-u)^2 进行高斯积分
    input: 四结点值u，i,j 单元格坐标，亦是单元格左下结点指标
    output: s:(uij-uxy)^2
    """
    p1 = (2 + 3 ** 0.5) / 6
    p2 = (2 - 3 ** 0.5) / 6
    p3 = 1. / 6.                         # 高斯点处四个基函数可能取值
    d1 = h * (1 - 1 / (3 ** 0.5)) / 2
    d2 = h * (1 + 1 / (3 ** 0.5)) / 2    # 高斯点在原单元的映射相对左下差值
    s1 = (u1 * p1 + u2 * p3 + u3 * p2 + u4 * p3 - uxy(xi + d1, yj + d1)) ** 2
    s2 = (u1 * p3 + u2 * p2 + u3 * p3 + u4 * p1 - uxy(xi + d2, yj + d1)) ** 2
    s3 = (u1 * p3 + u2 * p1 + u3 * p3 + u4 * p2 - uxy(xi + d2, yj + d2)) ** 2
    s4 = (u1 * p2 + u2 * p3 + u3 * p1 + u4 * p3 - uxy(xi + d1, yj + d2)) ** 2
    s = (s1 + s2 + s3 + s4) / 4
    return s


def uxy(x, y):
    """泊松方程真解"""
    return sin(pi * x) * sin(pi * y)


a_x = 0.                    # x方向边界
b_x = 1.
c_y = 0.                    # y方向边界
d_y = 1.
for m in range(4, 8):
    n = 2 ** m              # 剖分数
    h = (b_x - a_x) / n     # 步长

    A0 = coefficientmatrix(n)                # 计算系数矩阵
    f_1d = Matrix(zero_mat((n-1)*(n-1), 1))  # 代求向量，初始化
    b_matrix = fij(n, h, a_x, c_y)           # 计算右端矩阵b
    f_1d = cg(A0, b_matrix, f_1d)            # 共轭梯度法求解

    uij = []                                 # 转换成二维矩阵
    fT = f_1d.T().matrix[0]
    zero = [0 for i in range(n + 1)]
    uij.append(zero)
    for i in range(n - 1):
        uij.append([0] + fT[i * (n - 1):(i + 1) * (n - 1):] + [0])
    uij.append(zero)

    El2 = e_l2(uij, n, h, a_x, c_y)          # 计算误差
    print('n=', n, ' 误差： ', El2)

    plt.figure()
    x = [a_x + i * h for i in range(0, n + 1)]
    y = [c_y + i * h for i in range(0, n + 1)]
    plt.contourf(x, y, uij, n, cmap='seismic')
    plt.colorbar()
    plt.savefig(str(m - 3) + '.png')
plt.show()

