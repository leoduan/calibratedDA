from sympy import *
x = Symbol('x')
t = Symbol('t')
s = Symbol('s')

print integrate( exp(- 0.1*t**2)/t, (t, 0, 0.1))