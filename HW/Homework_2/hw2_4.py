import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['axes.unicode_minus'] = False


def f(x, r):
    return r * x * (1 - x)


def xr1(r):
    if r >= -1 and r <= 1:
        return 0
    elif r > 1 and r <= 3:
        return 1 - 1 / r


def xr2(r):
    # 取初值为0.7，最大迭代次数1000，迭代截断e=1e-4
    x = [0.7, 0, 0, 0]
    e = 1e-5
    for i in range(1001):
        x[1] = f(x[0], r)
        x[2] = f(x[1], r)
        x[3] = f(x[2], r)
        if i == 1000:
            print("r= ", r, "未收敛")
            return [False]
        elif abs(x[3] - x[1]) < e and abs(x[2] - x[0]) < e:
            return [min(x), max(x)]
        x[0] = f(x[3], r)


x1 = [-1 + 0.001 * i for i in range(0, 4001)]
y1 = [xr1(r) for r in x1]

x2 = [3]
x = xr2(x2[-1])
y3 = [x[0]]
y4 = [x[1]]
while True:
    x2.append(x2[-1] + 0.01)
    x = xr2(x2[-1])
    if isinstance(x[0], bool):
        x2.pop(-1)
        break
    y3.append(x[0])
    y4.append(x[1])
    if abs(x2[-1]**2 * (1 - 2 * y3[-1]) * (1 - 2 * y4[-1])) == 1:
        break
    elif abs(x2[-1]**2 * (1 - 2 * y3[-1]) * (1 - 2 * y4[-1])) > 1:
        y3.pop(-1), y4.pop(-1), x2.pop(-1)

xlabel1 = '$r$'
ylabel1 = '$x^*$'
plt.figure(figsize=(7.5, 5), dpi=100)  # 图片初始化与大小设置
plt.rcParams['font.sans-serif'] = ['SimHei']  # 中文显示设置
plt.plot(x1, y1, linestyle='-', marker='', linewidth=1)  # 绘图
plt.plot(x2, y3, linestyle='-', marker='', linewidth=1)
plt.plot(x2, y4, linestyle='-', marker='', linewidth=1)
plt.xlabel(xlabel1, size=12)  # 设置坐标轴标签
plt.ylabel(ylabel1, size=12)  # 设置坐标轴标签
plt.title('$x^*-r$', fontproperties="SimHei", size=15)  # 设置标题
plt.savefig('Homework_2/4.png')

plt.show()
