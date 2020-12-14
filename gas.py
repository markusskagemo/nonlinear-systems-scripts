from sympy import diff, simplify, expand, Q, ask, solve
from sympy import latex, symbols, exp, cos, sin, re, atan, pi
from sympy import Function, Symbol
from sympy import S
from sympy.parsing.latex import parse_latex
from sympy.matrices import Matrix

import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'


x_1, x_2, x_1dot, x_2dot, u, e = symbols('x_1 x_2 x_1dot x_2dot u e', real=True)
t = Symbol('t', positive=True, real=True)
V, V_dot = symbols('V V_dot', cls=Function)

    
def lasalle(x_1dot, x_2dot, V_dot):
    ls_x_1dot = x_1dot
    ls_x_2dot = x_2dot
    ls_sols = set()
    ls_unforced_sols = set()
    for x, s in solve(V_dot, (x_1, x_2)):
        if x == x_1 and x in x_1dot.free_symbols:
            ls_x_1dot = x_1dot.subs(x, s)
        elif x == x_2 and x in x_2dot.free_symbols:
            ls_x_2dot = x_2dot.subs(x, s)
        for sol in solve([ls_x_1dot, ls_x_2dot], (x_1, x_2)):
            ls_sols.add(str(sol))
        for sol in solve([ls_x_1dot.subs(u, 0), ls_x_2dot.subs(u, 0)], (x_1, x_2)):
            ls_unforced_sols.add(str(sol))
    return ls_sols, ls_unforced_sols


def eq_points(x_1dot, x_2dot):
    eqp = None
    unforced_eqp = None
    unforced_zero_eq = False
    
    free_1 = []
    if hasattr(x_1dot, 'free_symbols'):
        free_1 = x_1dot.free_symbols
    free_2 = []
    if hasattr(x_2dot, 'free_symbols'):
        free_2 = x_2dot.free_symbols
    
    if u in [*free_1, *free_2]:
        unforced_eqp = solve([x_1dot.subs(u, 0), x_2dot.subs(u, 0)], (x_1, x_2))
        if type(unforced_eqp) is dict:
            unforced_eqp = [tuple(list(unforced_eqp.values()))]
        print('Unforced system eq. points:', unforced_eqp, end='\n\n')
        if '(0, 0)' in [str(eq) for eq in unforced_eqp]:
            unforced_zero_eq = True
            print('Unforced eq in origin')
    else:
        eqp = solve([x_1dot, x_2dot], (x_1, x_2))
        if type(eqp) is dict:
            if x_1 not in eqp.keys():
                eqp[x_1] = S.Zero
            if x_2 not in eqp.keys():
                eqp[x_2] = S.Zero
            eqp = [tuple(list(eqp.values()))]
            
        print('System eq. points:', eqp)
    return eqp, unforced_eqp, unforced_zero_eq


def jacobian(x_1dot, x_2dot):
    x = Matrix([x_1, x_2])
    return Matrix([x_1dot, x_2dot]).jacobian(x)
    
    
def classify_eq_points(eqp, J):
    for eq_x_1, eq_x_2 in eqp:
        eig_mat = J.subs({x_1: eq_x_1, x_2: eq_x_2, u: 0})
        print(f'Linearized sys. for ({eq_x_1}, {eq_x_2}): A = {eig_mat}')
        eigs_d = eig_mat.eigenvals()
        eigs = eigs_d.keys()
        print(f'Eigenvals for ({eq_x_1, eq_x_2}): {list(eigs)}\nMultiplicity: {list(eigs_d.values())})')
    
        if any(not eig.is_real for eig in eigs):
            complex = True
            node = 'focus'
            saddle = 'Center'
        else:
            complex = False
            node = 'node'
            saddle = 'Saddle'
        
        if any(t in eig.free_symbols for eig in eigs):
            print('Time varying eigenvalues (t). Cannot determine truth value of Relational')
            print('Manually classify eigs:', eigs)
            continue
        try:
            eigs = sorted([re(eig) for eig in eigs])
        except Exception as e:
            print(f'Exception: {e}.\nSkipping eq. point ({eq_x_1}, {eq_x_2})')
            continue
            
        if len(eigs) == 1:
            if eigs[0] < 0:
                print(f'Class: Stable {node}.\n')
            elif eigs[0] > 0:
                print(f'Class: Unstable {node}.\n')
        elif len(eigs) == 2:
            if complex:
                if eigs[0] < 0:
                    print(f'Class: Stable {node}.\n')
                elif eigs[0] > 0:
                    print(f'Class: Unstable {node}.\n')
                else:
                    print(f'Class: {saddle} point.\n')
            else:
                print(eigs)
                if eigs[0] < eigs[1] < 0:
                    print(f'Class: Stable {node}.\n')
                elif eigs[1] > eigs[0] > 0:
                    print(f'Class: Unstable {node}.\n')
                else:
                    print(f'Class: {saddle} point.\n')

    
if __name__ == '__main__':
    # Specify system and LF
    x_1dot = -x_1
    x_2dot = (x_1*x_2 - 1)*x_2**3 + (x_1*x_2 - 1 + x_1**2)*x_2
    V = x_1**2/2 + x_2**2/2
    
    # Classify equilibrium points
    eqp, unforced_eqp, uze = eq_points(x_1dot, x_2dot)
    J = jacobian(x_1dot, x_2dot)
    print('Jacobian:', simplify(J), end='\n\n')
    if eqp is not None:
        classify_eq_points(eqp, J)
    if unforced_eqp is not None:
        print('Classifying unforced eq. points.')
        classify_eq_points(unforced_eqp, J)
    
    # \dot{V} algebra for given system
    plt_texts = ['$V = %s$' % latex(V)]
    V_dot = diff(V, x_1) * x_1dot + diff(V, x_2) * x_2dot
    plt_texts.append('$\dot{V} = %s$' % latex(V_dot))
    V_dot = simplify(V_dot)
    plt_texts.append('$\dot{V} = %s$' % latex(V_dot))
    
    # Weak definiteness tests
    if len(solve(V_dot == 0, (x_1, x_2))) > 1:
        raise ValueError('V_dot is not negative def.', V_dot)
    if len(solve(V == 0, (x_1, x_2))) > 1:
        raise ValueError('V is not positive def.', V)
    if eqp is not None and '(0, 0)' in [str(eq) for eq in eqp]:
        print('The origin is stable.', end='\n\n')

    try:
        # LaSalle (ls). See H18 problem 2 solution
        ls_sols, ls_unforced_sols = lasalle(x_1dot, x_2dot, V_dot)

        print(f'Unique LaSalle solutions: {ls_sols}')
        print(f'Unique unforced LaSalle solutions: {ls_unforced_sols}')
        if len(ls_sols) == 1 and ('(0.0, 0.0)' in ls_sols or '(0, 0)' in ls_sols):
            print("""The only solution that can stay identically in S is the trivial solution,
        and the origin is therefore GAS, assuming radially unbounded V, \dot{V} < 0 and V > 0 (by Corollary 4.2).""")
        if uze:
            if len(ls_unforced_sols) == 1 and ('(0, 0)' in ls_unforced_sols or '(0.0, 0.0)' in ls_sols):
                print(f'Unique unforced LaSalle solutions: {ls_unforced_sols}')
                print('The system is GAS for u(t) = 0. Assuming BIBO the system is ISS. (4.21 ISS hint)')
            else:
                print('Since the system is not 0-GAS, it cannot be input-to-state-stable (ISS).')
    
    ax = plt.axes([0, 0, 0.1, 0.1])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')
    plt.text(0.4, 0.4, '\n\n'.join(plt_texts))
    plt.show()