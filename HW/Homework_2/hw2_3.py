import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['axes.unicode_minus'] = False


def f(x, r):
    return r * x * (1 - x)


def opreate(x_0, r, n):
    x = [x_0]
    for i in range(n):
        x.append(f(x[-1], r))
    return x


a = []  # r = 3+0.1
for i in range(10):
    a.append(opreate(0.1*i, 3.1, 20))

xx = [i for i in range(21)]
xlabel1 = '$n$'
ylabel1 = '$x_n$'

plt.figure(figsize=(10, 5), dpi=100)           # 图片初始化与大小设置
plt.rcParams['font.sans-serif'] = ['SimHei']    # 中文显示设置
for i in range(10):      # 绘图
    plt.plot(xx, a[i], linestyle='-', marker='.', linewidth=1)
plt.xlabel(xlabel1, size=12)    # 设置坐标轴标签
plt.ylabel(ylabel1, size=12)    # 设置坐标轴标签
plt.title('$x_{n+1} = 3.1 x_n(1-x_n)$', fontproperties="SimHei", size=15)   # 设置标题
plt.savefig('Homework_2/3.1.png')

plt.show()
