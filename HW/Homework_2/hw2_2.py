import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['axes.unicode_minus'] = False


def xr(r):
    if r >= -1 and r <= 1:
        return 0
    elif r > 1 and r <= 3:
        return 1-1/r


x = [-1 + 0.001 * i for i in range(0, 4001)]
y = [xr(r) for r in x]
xlabel1 = '$r$'
ylabel1 = '$x^*$'

plt.figure(figsize=(7.5, 5), dpi=100)           # 图片初始化与大小设置
plt.rcParams['font.sans-serif'] = ['SimHei']    # 中文显示设置      
plt.plot(x, y, linestyle='-', marker='', linewidth=1)    # 绘图
plt.xlabel(xlabel1, size=12)    # 设置坐标轴标签
plt.ylabel(ylabel1, size=12)    # 设置坐标轴标签
plt.title('$x^*-r$', fontproperties="SimHei", size=15)   # 设置标题
plt.savefig('Homework_2/2.png')

plt.show()
