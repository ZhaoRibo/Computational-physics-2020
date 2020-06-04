import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['axes.unicode_minus'] = False


def f(x, r):
    return r * x * (1 - x)


def operate(x, r, n):
    for i in range(n):
        x = f(x,r)
    return x


for i in range(401):
    r = 0.01 * i
    x = operate(0.7, r, 500)
    for j in range(200):
        x = f(x, r)

