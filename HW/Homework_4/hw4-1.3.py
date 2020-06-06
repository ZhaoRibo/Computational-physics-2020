"""初始化采用一种近似均匀算法"""
from math import pi, sin, cos
import math
import copy
from wheels import Matrix, zero_mat


def initialize_x0(n):
    """初始化电子位置.

    cos(theta)和phi坐标在各自范围内随机取值.

    Args:
        n: 电子数.

    Returns:
        x0: 初始化的位置分布.
    """
    x0 = Matrix(zero_mat(2 * n, 1))
    i = 0
    for point in splot(n):
        if i >= n:
            break
        x0[i][0] = point[0]
        x0[i + n][0] = point[1]
        i += 1
    return x0


class Spherical(object):
    '''球坐标系'''
    def __init__(self, radial=1.0, polar=0.0, azimuthal=0.0):
        self.radial = radial
        self.polar = polar
        self.azimuthal = azimuthal

    def printCoordinate(self):
        return self.azimuthal, self.polar


def splot(N):
    s = Spherical()
    n = int(math.ceil(math.sqrt((int(N) - 2) / 4)))
    azimuthal = 0.5 * math.pi / n
    for a in range(-n, n + 1):
        s.polar = 0
        size = (n - abs(a)) * 4 or 1
        polar = 2 * math.pi / size
        for i in range(size):
            yield s.printCoordinate()
            s.polar += polar
        s.azimuthal += azimuthal


def rij(thetai, thetaj, phii, phij):
    """i,j两点间距（球坐标）.

    Args:
        theta_i, theta_j, phi_i, phi_j: i,j点坐标.

    Returns:
        r_ij: 两点间距.
    """
    r_ij = (sin((thetai - thetaj) / 2)**2 + sin(thetai) * sin(thetaj) * (sin(
        (phii - phij) / 2)**2))**0.5 * 2
    return r_ij


def u(x):
    """势能函数.

    Args:
        x: (theta,phi)^T.

    Returns:
        uu: 势能.
    """
    assert isinstance(x, Matrix), 'x不是Matrix'

    uu = 0
    n = int(x.row / 2)
    for i in range(n):
        for j in range(i):
            uu += 1 / rij(x[i][0], x[j][0], x[i + n][0], x[j + n][0])
    return uu


def g_u(x):
    """u的梯度(解析值)

    为尽量减少计算与误差，直接代入梯度解析式计算梯度.

    Args:
        x: (theta,phi)^T.

    Returns:
        g: 梯度矢量.
    """
    assert isinstance(x, Matrix), 'x不是Matrix'

    n = int(x.row / 2)
    g = Matrix(zero_mat(2 * n, 1))
    for i in range(n):
        a = 0
        b = 0
        for j in range(n):
            if j != i:
                a -= (sin(x[i][0] - x[j][0]) + cos(x[i][0]) * sin(x[j][0]) *
                      (1 - cos(x[i + n][0] - x[j + n][0]))) / (rij(
                          x[i][0], x[j][0], x[i + n][0], x[j + n][0])**3)
                b -= (sin(x[i][0]) * sin(x[j][0]) *
                      sin(x[i + n][0] - x[j + n][0])) / (rij(
                          x[i][0], x[j][0], x[i + n][0], x[j + n][0])**3)
        g[i][0] = a
        g[i + n][0] = b
    return g


def PRCG(x_0, maxk=1e4, ee=1e-12):
    """共轭梯度法:Polak-Ribiere Algorithm.

    Args:
        x_0: 初始点; maxk: 最大迭代次数; ee: 判停标准.

    Returns：
        x: 最优点; val: 最优解; k: 迭代次数.
    """
    x = copy.deepcopy(x_0)
    n = x.row  # 解维数
    k = 0      # 迭代次数
    g_0 = g_u(x)
    g0 = (g_0.T() * g_0) ** 0.5  # 初始梯度模长
    g = copy.deepcopy(g_0)       # x0处梯度初始化
    d = -g                       # 搜索方向初始化
    mode_g = g.T() * g           # 梯度模长初始化，减小计算量

    while k < maxk:
        if mode_g < ee:         # 判停
            break

        alpha = kiefer_gold(u, x, d, alpham=(1 / g0))   # 线搜索步长

        x = x + alpha * d           # x更新
        g_plus = g_u(x)             # g_(k+1)更新
        mode_gp = g_plus.T()*g_plus
        beta_PR = (mode_gp - g_plus.T() * g) / mode_g  # beta更新
        d = -g_plus + beta_PR * d   # 搜索方向更新
        g = copy.deepcopy(g_plus)   # 梯度更新
        mode_g = mode_gp
        k += 1

    for i in range(int(n / 2)):     # 0 <= theta <= pi, 0 <= phi < 2pi
        while True:
            if x[i][0] >= 0 and x[i][0] <= pi and x[
                    i + int(n / 2)][0] >= 0 and x[i + int(n / 2)][0] < 2 * pi:
                break
            if x[i][0] < 0:
                x[i][0] += pi
            elif x[i][0] > pi:
                x[i][0] -= pi
            if x[i + int(n / 2)][0] < 0:
                x[i + int(n / 2)][0] += 2 * pi
            elif x[i + int(n / 2)][0] >= 2 * pi:
                x[i + int(n / 2)][0] -= 2 * pi

    return x, u(x), k


def kiefer_gold(f, x, d, alpham=1., eps=1e-3):
    """线搜索kiefer黄金分割法.

    Args:
        f: 函数本体; x: 搜索初始点; d: 搜索方向; alpham: 搜索上限; eps: 判停标准.

    Returns:
        alpha: 搜索步长.
    """
    tau = (5 ** 0.5-1) / 2
    a = 0
    b = alpham
    x1 = a + (1 - tau) * (b - a)
    f1 = f(x + x1 * d)
    x2 = a + tau * (b - a)
    f2 = f(x + x2 * d)
    while True:
        if (b - a) < eps:
            if abs(b - alpham) <= eps:
                alpham = alpham * 2
                b = alpham
                x1 = a + (1 - tau) * (b - a)
                f1 = f(x + x1 * d)
                x2 = a + tau * (b - a)
                f(x + x2 * d)
            else:
                break
        if f1 > f2:
            a = x1
            x1 = x2
            f1 = f2
            x2 = a + tau * (b - a)
            f2 = f(x + x2 * d)
        else:
            b = x2
            x2 = x1
            f2 = f1
            x1 = a + (1 - tau) * (b - a)
            f1 = f(x + x1 * d)

    return (a + b) / 2


file = open('data_Q1.txt', 'w')
for n in range(2, 65):
    x0 = initialize_x0(n)
    xx, uu, kk = PRCG(x0, 500, 1e-4)
    k = 0
    while True:
        if kk < 500:
            k += 1
            if k == 5:
                break
        x0 = initialize_x0(n)
        xx1, uu1, kk1 = PRCG(x0, 500, 1e-4)
        if uu1 < uu:
            xx = copy.deepcopy(xx1)
            uu, kk = uu1, kk1
    xx, uu, kk2 = PRCG(xx, 1e4, 1e-12)
    kk = kk + kk2
    file.write('电子数: ' + str(n) + '\n')
    file.write('分布:\n' + str(xx.matrix) + '\n')
    file.write('势能: ' + str(uu) + '\n')
    file.write('迭代次数: ' + str(kk) + '\n\n')
    file.flush()
    print('电子数：', n)
    print('势能：', uu)
    print('迭代次数：', kk)
file.close()
