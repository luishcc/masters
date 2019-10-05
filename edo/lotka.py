import scipy as sp
import matplotlib.pyplot as plt

dt = 0.0001
time = 1000000

x = sp.zeros(time)
y = sp.zeros(time)

alpha = 0.6666
beta = 1.3333
delta = 1
gamma = 1

x[0] = 0.9
y[0] = 0.9

for t in range(1, time):
	print t, "/", time
	x[t] = (alpha*x[t-1] - beta * x[t-1] * y[t-1]) * dt + x[t-1]
	y[t] = (delta * x[t-1] * y[t-1] - gamma * y[t-1]) * dt + y[t-1]


plt.figure(1)
plt.plot(x, y, 'b.')

plt.figure(2)
plt.plot(sp.linspace(0,time*dt, time), x, 'k-')
plt.plot(sp.linspace(0,time*dt, time), y, 'b-')
plt.show()

