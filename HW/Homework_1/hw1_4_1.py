''' 三角形交流网络 '''
from wheels import Matrix
from wheels import zero_mat
from wheels import up_tri_msolve
import cmath
import time
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['axes.unicode_minus'] = False


def get_G(n, w):
    '''
    获取三角形交流网络的导纳矩阵
    取边上节点数为n，节点标号方式为从右下一个顶点开始，算作第一层，一次向上取第二层，第三层……
    return Matrix
    '''
    g = Matrix(zero_mat(int((n ** 2 + n) / 2), int((n ** 2 + n) / 2)))
    g[0][0] = 1j * (w - 1 / w)
    g[0][1] = -1j * w
    g[0][2] = 1j / w
    i = 1
    for k in range(1, n-1):
        g[i][i - k] = -1j * w
        g[i][i + 1] = -1
        g[i][i + k + 1] = -1j * w
        g[i][i + k + 2] = 1j / w
        g[i][i] = 1 + 1j * (2 * w - 1 / w)
        for j in range(1, k):
            s = i + j
            g[s][s - k - 1] = 1j / w
            g[s][s - k] = -1j * w
            g[s][s - 1] = -1
            g[s][s + 1] = -1
            g[s][s + k + 1] = -1j * w
            g[s][s + k + 2] = 1j / w
            g[s][s] = 2 + 2j * (w - 1 / w)
        s = i + k
        g[s][s - k - 1] = 1j / w
        g[s][s - 1] = -1
        g[s][s + k + 1] = -1j * w
        g[s][s + k + 2] = 1j / w
        g[s][s] = 1 + 1j * (w - 2 / w)
        i = s + 1
    g[i][i - n + 1] = -1j * w
    g[i][i + 1] = -1
    g[i][i] = 1 + 1j * w
    for j in range(1, n - 1):
        s = i + j
        g[s][s - n] = 1j / w
        g[s][s - n + 1] = -1j * w
        g[s][s - 1] = -1
        g[s][s + 1] = -1
        g[s][s] = 2 + 1j * (w - 1 / w)
    s = i + n - 1
    g[s][s - n] = 1j / w
    g[s][s - 1] = -1
    g[s][s] = 1 - 1j / w
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

    g = g.T()
    # 上三角矩阵回代得解
   
    # 直接调用wheels通法
    x = up_tri_msolve(g, I)
    return x[b-1][0]


start = time.time()
f = [0.001 * i for i in range(1, 1000)]+[0.001 * i for i in range(1001, 2000)]
R = []
arg = []
for w in f:
    g = get_G(5, w)
    x = solve_R(g, 10, 14, 5)
    R.append(abs(x))
    arg.append(cmath.phase(x))
end = time.time()
print('datatime: ', end - start, '\n')

xlabel = '$\omega = 2\pi f$'
ylabel1 = '$R$'
ylabel2 = '$\\varphi$'

plt.figure(figsize=(7.5, 5), dpi=100)
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.plot(f, R, linestyle='-', marker='', linewidth=1)
plt.xlabel(xlabel, size=12)         # 设置坐标轴标签
plt.ylabel(ylabel1, size=12)
plt.title('阻抗大小频率响应曲线(直接法）', fontproperties="SimHei", size=15)   # 设置标题
plt.savefig('4.1.1.png')

plt.figure(figsize=(7.5, 5), dpi=100)
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.plot(f, arg, linestyle='-', marker='', linewidth=1)
plt.xlabel(xlabel, size=12)         # 设置坐标轴标签
plt.ylabel(ylabel2, size=12)
plt.title('阻抗相角频率响应曲线(直接法)', fontproperties="SimHei", size=15)   # 设置标题
plt.savefig('4.1.2.png')

plt.show()
