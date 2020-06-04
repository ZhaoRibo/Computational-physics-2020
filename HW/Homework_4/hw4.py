from math import pi, sin, cos
import copy
import random
from wheels import Matrix, zero_mat, identity_mat, up_tri_msolve, low_tri_msolve


def initialize_x0(n):
    """随机初始化x0，使得初始各电子分开，避免出现分母为零的情况
    Args:
        电子数 n
    Returns:
        Matrix x0
    """
    x0 = Matrix(zero_mat(2 * n, 1))
    for i in range(n):
        x0[i][0] = random.uniform(0, pi)
        x0[i + n][0] = random.uniform(0, 2 * pi)
    return x0


def rij(thetai, thetaj, phii, phij):
    """i,j两点间距（球坐标）
    Args:
        theta_i, theta_j, phi_i, phi_j
    """
    r_ij = (sin((thetai - thetaj) / 2)**2 + sin(thetai) * sin(thetaj) * (sin(
        (phii - phij) / 2)**2))**0.5 * 2
    return r_ij


def u(x):
    """势能函数
    Args:
        x=(theta,phi)^T
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
    Args:
        x
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


def solve(A, b):
    """Ax = b 求解, 改进算法的Cholesky分解"""
    tl = copy.deepcopy(A)
    n = A.row
    for i in range(n):
        for j in range(i):
            for k in range(j):
                tl[i][j] -= tl[i][k] * tl[k][j]
            tl[j][i] = tl[i][j] / tl[j][j]
        for k in range(i):
            tl[i][i] -= tl[i][k] * tl[k][i]
    y = low_tri_msolve(tl, b)
    for i in range(n):
        tl[i][i] = 1
    x = up_tri_msolve(tl, y)
    return x


def PRCG(x_0, maxk=1e5, ee=1e-14):
    """共轭梯度法:Polak-Ribiere Algorithm
    Args:
        x_0 初始点，maxk 最大迭代次数，ee 判停标准
    Returns：
        x 最优点，val 最优解，k 迭代次数
    """
    x = copy.deepcopy(x_0)
    n = x.row  # 解维数
    k = 0       # 迭代次数
    kk = 0      # 重置监视
    g_0 = g_u(x) 
    g0 = (g_0.T() * g_0)** 0.5  # 初始梯度模长
    g = copy.deepcopy(g_0)  # x0处梯度初始化
    d = -g      # 搜索方向初始化

    while k < maxk:
        if (g.T() * g) < ee: break  # 判停
        if kk == n:                 # 重启判定
            d = -g
            kk = 0
        alpha = kiefer_gold(u, x, d, alpham=(1 / g0))   # 线搜索步长
        x = x + alpha * d
        g_plus = g_u(x)
        beta_PR = (g_plus.T() * (g_plus - g)) / (g.T() * g)
        d = -g_plus + beta_PR * d   # 搜索方向更新
        g = g_plus                  # 梯度更新
        kk += 1
        k += 1

    for i in range(int(n / 2)):  # 0 <= theta <= pi, 0 <= phi < 2pi
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


def kiefer_gold(f, x, d, alpham=1, eps=1e-10):
    """线搜索kiefer黄金分割法
    Args:
        f 函数本体，x 搜索初始点，d 搜索方向，alpham 搜索上限，eps 判停标准
    Returns:
        alpha 步长
    """
    tau = (5 ** 0.5-1) / 2
    a = 0
    b = alpham
    x1 = a + (1 - tau) * (b - a); f1 = f(x + x1 * d)
    x2 = a + tau * (b - a); f2 = f(x + x2 * d)
    while True:
        if (b - a) < eps:
            if abs(b - alpham) <= eps:
                alpham = alpham * 2
                b = alpham
                x1 = a + (1 - tau) * (b - a); f1 = f(x + x1 * d)
                x2 = a + tau * (b - a); f(x + x2 * d)
            else:
                break
        if f1 > f2:
            a = x1; x1 = x2; f1 = f2
            x2 = a + tau * (b - a); f2 = f(x + x2 * d)
        else:
            b = x2; x2 = x1; f2 = f1;
            x1 = a + (1 - tau) * (b - a); f1 = f(x + x1 * d)
    
    return (a + b) / 2


for n in range(2, 65):
    x0 = initialize_x0(n)
    xx, uu, kk = PRCG(x0)
    for i in range(4):
        x0 = initialize_x0(n)
        xx1, uu1, kk1 = PRCG(x0)
        if uu1 < uu:
            xx, uu, kk = xx1, uu1, kk1
    print(xx.matrix)
    print('势能：', uu)
    print('迭代次数：', kk)
