

import sympy

x = sympy.Symbol('x')
y = sympy.Symbol('y')
z = sympy.Symbol('z')


norm3 = sympy.sqrt( x**2 + y**2 + z**2 )

print 'norm3 = ',norm3
print 'd(norm3)/dx = ',sympy.diff( norm3, x )

a0 = sympy.Symbol('a0')
a1 = sympy.Symbol('a1')
a2 = sympy.Symbol('a2')
b0 = sympy.Symbol('b0')
b1 = sympy.Symbol('b1')
b2 = sympy.Symbol('b2')

dot3 = a0*b0 + a1*b1 + a2*b2

print 'dot3 =',dot3
print 'd(dot3)/da0 =', sympy.diff( dot3, a0 )



