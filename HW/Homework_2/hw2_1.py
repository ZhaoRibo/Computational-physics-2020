import matplotlib.pyplot as plt


def f(x, r):
    return r * x * (1 - x)


def opreate(x_0, r, n):
    x = [x_0]
    for i in range(n):
        x.append(f(x[-1], r))
    return x


a = []  # r = 0.5
b = []  # r = 1.5
for i in range(11):
    a.append(opreate(0.1 * i, 0.5, 10))
    b.append(opreate(0.1 * i, 1.5, 10))

xx = [i for i in range(11)]
xlabel1 = '$n$'
ylabel1 = '$x_n$'

plt.figure(figsize=(7.5, 5), dpi=100)           # 图片初始化与大小设置
plt.rcParams['font.sans-serif'] = ['SimHei']    # 中文显示设置
for i in range(11):      # 绘图
    plt.plot(xx, a[i], linestyle='-', marker='.', linewidth=1)
plt.xlim([0, 11])   # 设置刻度
plt.ylim([0, 1])    # 设置刻度
plt.xlabel(xlabel1, size=12)    # 设置坐标轴标签
plt.ylabel(ylabel1, size=12)    # 设置坐标轴标签
plt.title('$x_{n+1} = 0.5 x_n(1-x_n)$', fontproperties="SimHei", size=15)   # 设置标题
plt.savefig('Homework_2/1.1.png')

plt.figure(figsize=(7.5, 5), dpi=100)           # 图片初始化与大小设置
plt.rcParams['font.sans-serif'] = ['SimHei']    # 中文显示设置
for i in range(11):      # 绘图
    plt.plot(xx, b[i], linestyle='-', marker='.', linewidth=1)
plt.xlim([0, 11])   # 设置刻度
plt.ylim([0, 1])    # 设置刻度
plt.xlabel(xlabel1, size=12)    # 设置坐标轴标签
plt.ylabel(ylabel1, size=12)    # 设置坐标轴标签
plt.title('$x_{n+1} = 1.5x_n(1-x_n)$', fontproperties="SimHei", size=15)   # 设置标题
plt.savefig('Homework_2/1.2.png')

plt.show()
