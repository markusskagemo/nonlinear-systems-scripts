## nonlinear-systems-scripts
Verify inferred function components:  
f = -x_1**3 + exp(x_1) - 1, g = 1, h = -x_2**3  

```math
\phi(x_1) = k_{1} \\left(x_{1}^{4} - x_{1} e^{x_{1}} + x_{1}\\right) \newline
\dot{\\phi}(x_1) = k_{1} \\left(4 x_{1}^{3} - x_{1} e^{x_{1}} - e^{x_{1}} + 1\\right) \newline
  
u(x_1, x_2) = - k_{1} \\left(x_{1}^{3} - x_{2} - e^{x_{1}} + 1\\right) \\left(4 x_{1}^{3} - x_{1} e^{x_{1}} - e^{x_{1}} + 1\\right) - k_{2} z - x_{1} + x_{2}^{3} \newline
u(x_1, x_2) = - k_{1} \\left(x_{1}^{3} - x_{2} - e^{x_{1}} + 1\\right) \\left(4 x_{1}^{3} - x_{1} e^{x_{1}} - e^{x_{1}} + 1\\right) + k_{2} \\left(k_{1} x_{1} \\left(x_{1}^{3} - e^{x_{1}} + 1\\right) - x_{2}\\right) - x_{1} + x_{2}^{3} \newline
  
V_2 = \\frac{x_{1}^{2}}{2} + \\frac{z^{2}}{2} \newline
V_2 = \\frac{x_{1}^{2}}{2} + \\frac{\\left(- k_{1} \\left(x_{1}^{4} - x_{1} e^{x_{1}} + x_{1}\\right) + x_{2}\\right)^{2}}{2} \newline
\dot{V}_2 = - k_{2} z^{2} + x_{1} \\left(k_{1} x_{1} \\left(x_{1}^{3} - e^{x_{1}} + 1\\right) - x_{1}^{3} + e^{x_{1}} - 1\\right)
```
