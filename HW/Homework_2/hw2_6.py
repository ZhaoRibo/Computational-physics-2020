import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['axes.unicode_minus'] = False


def f(x, r):
    return r * x * (1 - x)


def operate(x, r, n):
    for i in range(n):
        x = f(x, r)
    return x


plt.figure(figsize=(10, 5), dpi=100)           # 图片初始化与大小设置
plt.rcParams['font.sans-serif'] = ['SimHei']    # 中文显示设置

for i in range(401):
    r = 0.01 * i
    x = operate(0.7, r, 300)
    for j in range(100):
        x = f(x, r)
        plt.scatter(r, x, s=1)

# plt.plot(xx, a[1], linestyle='-', marker='', linewidth=1)
plt.xlabel('$r$', size=12)    # 设置坐标轴标签
plt.ylabel('$x^*$', size=12)    # 设置坐标轴标签
plt.title('$x^*-r$', fontproperties="SimHei", size=15)   # 设置标题
plt.savefig('Homework_2/6.1.png')


plt.show()
