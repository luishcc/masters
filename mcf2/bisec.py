import math as mt

p=0.4

def f(g):
	return (p+g)*(1+(3*((g-p)**2/(g+p)**2))/(10+ mt.sqrt(4-3*((g-p)**2/(g+p)**2))))-4*p*g
	
def bisection(a,b,tol):
	c = (a+b)/2.0
	while (b-a)/2.0 > tol:
		if f(c) == 0:
			return c
		elif f(a)*f(c) < 0:
			b = c
		else :
			a = c
		c = (a+b)/2.0
		
	return c


print bisection(0, 3, 0.0000001)
