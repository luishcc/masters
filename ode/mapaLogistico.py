import scipy as sp
import matplotlib.pyplot as plt

mi = 0.950     # 0 < Mi <1
x_0 = 0.25    # 0 < x_0 <1
x2_0 = 0.2501    # 0 < x_0 <1



tempo = 40
n = sp.linspace(0, tempo, tempo)


x = sp.zeros(tempo)
x2 = sp.zeros(tempo)

x_para = sp.linspace(0, 1, 1000)
x_ret = sp.linspace(0, 1, 1000)
y_para = sp.zeros(1000)
for i in range(1000):
    y_para[i] = 4*mi*x_para[i]*(1-x_para[i])



x[0] = x_0
x2[0] = x2_0
for i in range(1, tempo):
    x[i] = 4*mi*x[i-1]*(1-x[i-1])
    x2[i] = 4 * mi * x2[i - 1] * (1 - x2[i - 1])
    plt.figure(1)
    plt.plot(x[i-1], x[i], 'bo')
    plt.plot(x2[i-1], x2[i], 'ko')




plt.figure(2)
plt.xlim(0, tempo)
plt.ylim(0, 1)
plt.plot(n, x, 'b-')
plt.plot(n, x2, 'k-')


plt.figure(1)
plt.plot(x_para, y_para, 'k-')
plt.plot(x_ret, x_ret, 'k-')
plt.xlim(0,1)
plt.ylim(0,1)

plt.show()


