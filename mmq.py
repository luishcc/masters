import numpy as np
import matplotlib.pyplot as plt

n = 10
a = 5
x = np.linspace(0, a, n)
y = np.zeros(n)
for i in xrange(n):
	y[i] = np.random.normal(x[i],1)

def estCoef(_x, _y, _n):
	mean_x, mean_y = np.mean(_x),  np.mean(_y)

	SS_xy = np.sum(_y * _x - _n * mean_y * mean_x)
	SS_xx = np.sum(_x * _x - _n * mean_x * mean_x)

	b_1 = SS_xy / SS_xx
	b_0 = mean_y - b_1 * mean_x
	
	return(b_0, b_1)

b = estCoef(x, y, n)
y_pred = b[0] + b[1] * x

plt.figure(1)
plt.plot(x, y, 'bo')
plt.plot(x, y_pred, 'k-')
plt.show()



