from sympy.parsing.latex import parse_latex
from sympy import symbols, diff, exp, simplify, expand, latex, solve, atan
from sympy import Function, Wild

import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 6})

u, z, x_1, x_2, z = symbols('u z x_1 x_2, z', real=True)
k_1, k_2 = symbols('k_1 k_2', positive=True, real=True)
f, g, h, V_1, V_2, phi = symbols('f g h V_1 V_2 phi')


if __name__ == '__main__':
    # Specify system and LF
    x_1dot = x_1/3 + x_2
    x_2dot = 0#-x_2**2 + u
    V_1 = x_1**2/2
    
    # Infer function components f, g, h from system
    g_unsub = x_1dot.as_independent(x_2)[1]
    g = g_unsub.subs({x_2: 1})
    f = x_1dot - g_unsub
    h = x_2dot - x_2dot.as_independent(u)[1]
    print('Verify inferred function components:')
    print('f = {}\ng = {}\nh = {}'.format(f, g, h))
    
    V_1dot = diff(V_1, x_1) * x_1dot
    # Choose phi so \dot{V_1} < 0
    phi = (k_1 + 1) * solve(V_1dot, x_2)[0]
    print('Verify inferred \phi:', phi)
    phi_dot = diff(phi, x_1)
    
    
    # Verify that phi holds conditions
    assert x_2 not in phi.free_symbols, 'x_2 should not be in \phi(x_1)'
    assert x_2 not in phi_dot.free_symbols, 'x_2 should not be in \dot{\phi}(x_1)'
    assert phi.subs({x_1: 0, x_2: 0}) == 0, f'phi(0) != 0. phi(0) = {phi.subs({x_1: 0, x_2: 0})}'
    
    inp = -k_2 * z - (x_1 * g + h - phi_dot * (f + g * z)) # Controlling input
    #z_dot = h + u - phi_dot * x_1dot * (f + g * z)
    new_x_1dot = simplify(f + g * phi + g * z)
    z_dot = x_2dot - phi_dot * new_x_1dot #x_1dot.replace(x_2 - phi, z)
    #z_dot = diff(x_2 - phi, x_1) * x_1dot + diff(x_2 - phi, x_2) * x_2dot
    V_2 = V_1 + V_1.subs(x_1, z)
    V_2dot = x_1 * (f + g * phi) - k_2 * z**2
    #u = u.subs(z, x_2 - phi)
    V_2dot = V_2dot.subs(u, inp)
    
    # Weak definiteness tests
    if len(solve(V_2dot == 0, (x_1, x_2, z))) > 1:
        raise ValueError('V_dot is not negative def.', V_dot)
    if len(solve(V_2 == 0, (x_1, x_2))) > 1:
        raise ValueError('V is not positive def.', V)
    
    plt_texts = ['$f(x_1) = {}$\n$g(x_1) = {}$\n$h(x_2) = {}$'.format(latex(f), latex(g), latex(h))]
    plt_texts.append('$\phi(x_1) = %s$' % latex(phi))
    plt_texts.append('$\dot{\phi}(x_1) = %s$\n' % latex(phi_dot))
    plt_texts.append('Substituting $z = x_2 - \phi(x_1)$, which gives the new system')
    plt_texts.append('$\dot{x_1} = f(x_1) + g(x_1)\phi(x_1) + g(x_1)z$')
    plt_texts.append('$\dot{x_1} = %s$' % latex(new_x_1dot))
    plt_texts.append('$\dot{z} = h(x_2) + u - \partial\phi(x_1)/\partial x_1(f(x_1) +g(x_1)z)$')
    plt_texts.append('$\dot{z} = %s$\n' % latex(simplify(z_dot)))
    plt_texts.append('Choosing the control input $u$')
    plt_texts.append('$u = %s$' % latex(simplify(inp)))
    plt_texts.append('$u = %s$\n' % latex(simplify(inp.subs(z, x_2 - phi))))
    plt_texts.append('$V_2 = %s$' % latex(simplify(V_2)))
    plt_texts.append('$V_2 = %s$' % latex(simplify(V_2).subs(z, x_2 - phi)))
    plt_texts.append('$\dot{V}_2 = %s < 0$' % latex(simplify(V_2dot)))
    plt_texts.append("""
Since $V_2$ is positive definite and radially unbounded, this makes the origin $(x_1, z) = (0, 0)$ GAS.
Since $\phi(0) = 0$, the origin of the original system in $(x_1, x_2)$ coincides
with that of the transformed coordinates and is also GAS.
    """)
    
    # Render LaTeX
    ax = plt.axes([0, 0, 0.2, 0.2])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')
    plt.text(0.1, 0.1, '\n\n'.join(plt_texts))
    plt.show()