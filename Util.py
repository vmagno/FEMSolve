
from sympy import Symbol

x = Symbol('x')

def SymbolAsFunc(symbol):
    def _function(a):
        return symbol.subs(x, a)
    return _function
