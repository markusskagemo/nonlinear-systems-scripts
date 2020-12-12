from sympy.parsing.latex import parse_latex
from sympy import symbols, diff, exp, simplify, expand, latex, solve, atan
from sympy import Function, Wild

import matplotlib.pyplot as plt


u, z, x_1, x_2, z = symbols('u z x_1 x_2, z', real=True)
k_1, k_2 = symbols('k_1 k_2', positive=True, real=True)
f, g, h, V_1, V_2, phi = symbols('f g h V_1 V_2 phi')


if __name__ == '__main__':
    # Specify system and LF
    x_1dot = 1/exp(x_1) - 2*x_1**3 - 1 + x_2/x_1**2 + x_2*x_1**3
    x_2dot = -x_2**5 + u + 1 + exp(atan(x_2))
    V_1 = x_1**2/2
    
    # Infer function components f, g, h from system
    g_unsub = x_1dot.as_independent(x_2)[1]
    g = g_unsub.subs({x_2: 1})
    f = x_1dot - g_unsub
    h = x_2dot - x_2dot.as_independent(u)[1]
    print('Verify inferred function components:')
    print('f = {}\ng = {}\nh = {}'.format(f, g, h))
    
    V_1dot = diff(V_1, x_1) * (f + phi * g)
    # Choose phi so \dot{V_1} < 0
    phi = k_1 * solve(V_1dot < 0, phi).rhs # > TODO: Improve phi form guess
    phi_dot = diff(phi, x_1)
    # Verify that phi holds conditions
    assert x_2 not in phi.free_symbols, 'x_2 should not be in \phi(x_1)'
    assert x_2 not in phi_dot.free_symbols, 'x_2 should not be in \dot{\phi}(x_1)'
    assert phi.subs({x_1: 0, x_2: 0}) == 0, f'phi(0) != 0. phi(0) = {phi.subs({x_1: 0, x_2: 0})}'

    # Prepare substitution
    #z_dot = h - phi_dot * x_1dot #(f + g * x_2)
    
    # > Verify stabilizing input format
    inp = -k_2 * z - (x_1 * g + h - phi_dot * (f + g * x_2))

    V_2 = V_1 + V_1.subs(x_1, z)
    V_2dot = x_1 * (f + g * phi) - k_2 * z**2
    V_2dot = V_2dot.subs(u, inp)
    
    # Weak definiteness tests
    if len(solve(V_2dot == 0, (x_1, x_2, z))) > 1:
        raise ValueError('V_dot is not negative def.', V_dot)
    if len(solve(V_2 == 0, (x_1, x_2))) > 1:
        raise ValueError('V is not positive def.', V)

    plt_texts = ['$\phi(x_1) = %s$' % latex(phi)]
    plt_texts.append('$\dot{\phi}(x_1) = %s$' % latex(phi_dot))
    plt_texts.append('$u(x_1, x_2) = %s$' % latex(simplify(inp)))
    plt_texts.append('$u(x_1, x_2) = %s$' % latex(simplify(inp.subs(z, x_2 - phi))))
    plt_texts.append('$V_2 = %s$' % latex(simplify(V_2)))
    plt_texts.append('$V_2 = %s$' % latex(simplify(V_2).subs(z, x_2 - phi)))
    plt_texts.append('$\dot{V}_2 = %s$' % latex(simplify(V_2dot)))
    
    # Render LaTeX
    ax = plt.axes([0, 0, 0.1, 0.1])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')
    plt.text(0.4, 0.4, '\n\n'.join(plt_texts))
    plt.show()