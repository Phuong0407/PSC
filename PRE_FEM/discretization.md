For a function $\psi(x,y)$ define as $u(x,y)$ on a square grid, and $u(x, \beta x + y)$ on oblique grid:

===================================================================
===================================================================
For rectangular grid + oblique grid interface:

\[
\begin{aligned}
& u(x + \Delta x_2, y) - u(x, y) = \left( \dfrac{\partial u}{\partial x} + \beta\dfrac{\partial u}{\partial y} \right) \Delta x_2 \Longrightarrow \dfrac{u(x + \Delta x_2, y) - u(x, y)}{\Delta x_2} = \dfrac{\partial u}{\partial x} + \beta\dfrac{\partial u}{\partial y} \\
& \dfrac{u(x + \Delta x_2, y) - u(x, y)}{\Delta x_2} = \dfrac{\partial u}{\partial x} + \beta\dfrac{\partial u}{\partial y} + \dfrac{\Delta x_2}{2} \left( \dfrac{\partial^2 u}{\partial x^2} + 2\beta\dfrac{\partial^2 u}{\partial x \partial y} + \beta\dfrac{\partial^2 u}{\partial y^2} \right) \\
& \dfrac{u(x + \Delta x_1, y) - u(x, y)}{\Delta x_1} = - \dfrac{\partial u}{\partial x} + \dfrac{\Delta x_1}{2} \dfrac{\partial^2 u}{\partial x^2} \\
& \dfrac{u(x, y + \Delta y_1) - u(x, y)}{\Delta y_1} = \dfrac{\partial u}{\partial y} \\
& \dfrac{u(x, y + \Delta y_1) + u(x, y - \Delta y_1) - 2 u(x, y)}{\Delta y_1^2} = \dfrac{\partial^2 u}{\partial y^2}
\end{aligned}
\]

Therefore 1
\[
\begin{aligned}
\dfrac{u(x + \Delta x_2, y) - u(x, y)}{\Delta x_2} + \dfrac{u(x - \Delta x_1, y) - u(x, y)}{\Delta x_1} = \beta\dfrac{\partial u}{\partial y} + \dfrac{\Delta x_1 + \Delta x_2}{2} \dfrac{\partial^2 u}{\partial x^2} + \dfrac{\Delta x_2}{2} \beta \left( 2\dfrac{\partial^2 u}{\partial x \partial y} + \beta\dfrac{\partial^2 u}{\partial y^2} \right)
\end{aligned}
\]

We could see that

\[
\begin{aligned}
\dfrac{\partial^2 u}{\partial x \partial y} + \beta\dfrac{\partial^2 u}{\partial y^2} = \dfrac{\partial}{\partial y}\left( \dfrac{\partial u}{\partial x} + \beta\dfrac{\partial u}{\partial y} \right) = & \dfrac{1}{\Delta y_1} \left[ \left(\dfrac{\partial u}{\partial x} + \beta\dfrac{\partial u}{\partial y} \right)(x, y + \Delta y_1) - \left(\dfrac{\partial u}{\partial x} + \beta\dfrac{\partial u}{\partial y} \right)(x, y) \right] \\
= & \dfrac{u(x + \Delta x_2, y + \Delta y_1) - u(x, y + \Delta y_1)}{\Delta x_2 \Delta y_1} - \dfrac{u(x + \Delta x_2, y) - u(x, y)}{\Delta x_2 \Delta y_1}
\end{aligned}
\]

Therefore 2
\[
\begin{aligned}
\dfrac{u(x + \Delta x_2, y) - u(x, y)}{\Delta x_2} + \dfrac{u(x - \Delta x_1, y) - u(x, y)}{\Delta x_1} = & \beta\dfrac{\partial u}{\partial y} + \dfrac{\Delta x_2 + \Delta x_3}{2} \dfrac{\partial^2 u}{\partial x^2} + \dfrac{\Delta x_3}{2} \beta \left( 2\dfrac{\partial^2 u}{\partial x \partial y} + 2 \beta\dfrac{\partial^2 u}{\partial y^2} \right) - \dfrac{\Delta x_2}{2} \beta^2 \dfrac{\partial^2 u}{\partial y^2} \\
= & \beta\dfrac{u(x, y + \Delta y_1) - u(x, y)}{\Delta y_1} + \dfrac{\Delta x_1 + \Delta x_2}{2} \dfrac{\partial^2 u}{\partial x^2} \\
& + \Delta x_2 \beta \left( \dfrac{u(x + \Delta x_2, y + \Delta y_1) - u(x, y + \Delta y_1)}{\Delta x_2 \Delta y_1} - \dfrac{u(x + \Delta x_2, y) - u(x, y)}{\Delta x_2 \Delta y_1} \right) - \dfrac{\Delta x_2}{2} \beta^2 \dfrac{\partial^2 u}{\partial y^2} \\
\end{aligned}
\]

Therefore 3
\[
\begin{aligned}
2 \dfrac{u(x + \Delta x_2, y) - u(x, y)}{\Delta x_2 (\Delta x_1 + \Delta x_2)} + 2 \dfrac{u(x - \Delta x_1, y) - u(x, y)}{\Delta x_1 (\Delta x_1 + \Delta x_2)} - 2\beta\dfrac{u(x, y + \Delta y_1) - u(x, y)}{\Delta y_1 (\Delta x_1 + \Delta x_2)} - 2 \beta \left( \dfrac{u(x + \Delta x_2, y + \Delta y_1) - u(x, y + \Delta y_1)}{\Delta y_1 (\Delta x_1 + \Delta x_2)} - \dfrac{u(x + \Delta x_2, y) - u(x, y)}{\Delta y_1 (\Delta x_1 + \Delta x_2)} \right) = & \dfrac{\partial^2 u}{\partial x^2} - \dfrac{\Delta x_2}{\Delta x_1 + \Delta x_2} \beta^2 \dfrac{\partial^2 u}{\partial y^2}
\end{aligned}
\]

Therefore 4
\[
\begin{aligned}
& 2 \dfrac{u(x + \Delta x_2, y) - u(x, y)}{\Delta x_2 (\Delta x_1 + \Delta x_2)} + 2 \dfrac{u(x - \Delta x_1, y) - u(x, y)}{\Delta x_1 (\Delta x_1 + \Delta x_2)} - 2\beta\dfrac{u(x, y + \Delta y_1) - u(x, y)}{\Delta y_1 (\Delta x_1 + \Delta x_2)} - 2 \beta \left( \dfrac{u(x + \Delta x_2, y + \Delta y_1) - u(x, y + \Delta y_1)}{\Delta y_1 (\Delta x_1 + \Delta x_2)} - \dfrac{u(x + \Delta x_2, y) - u(x, y)}{\Delta y_1 (\Delta x_1 + \Delta x_2)} \right) \\
& + \left( \dfrac{\Delta x_2}{\Delta x_1 + \Delta x_2} \beta^2 + 1 \right) \dfrac{u(x, y + \Delta y_1) + u(x, y - \Delta y_1) - 2 u(x, y)}{\Delta y_1^2} = \dfrac{\partial^2 u}{\partial x^2} + \dfrac{\partial^2 u}{\partial y^2} = \Delta u
\end{aligned}
\]


===================================================================
===================================================================

For oblique grid + rectangular grid interface:
\[
\begin{aligned}
& u(x - \Delta x_2, y) - u(x, y) = - \left( \dfrac{\partial u}{\partial x} + \beta\dfrac{\partial u}{\partial y} \right) \Delta x_2 \Longrightarrow \dfrac{u(x - \Delta x_2, y) - u(x, y)}{\Delta x_2} = - \left( \dfrac{\partial u}{\partial x} + \beta\dfrac{\partial u}{\partial y} \right) \\
& \dfrac{u(x - \Delta x_2, y) - u(x, y)}{\Delta x_2} = - \left( \dfrac{\partial u}{\partial x} + \beta\dfrac{\partial u}{\partial y} \right) + \dfrac{\Delta x_2}{2} \left( \dfrac{\partial^2 u}{\partial x^2} + 2\beta\dfrac{\partial^2 u}{\partial x \partial y} + \beta\dfrac{\partial^2 u}{\partial y^2} \right) \\
& \dfrac{u(x + \Delta x_3, y) - u(x, y)}{\Delta x_3} = \dfrac{\partial u}{\partial x} + \dfrac{\Delta x_3}{2} \dfrac{\partial^2 u}{\partial x^2} \\
& \dfrac{u(x, y + \Delta y_3) - u(x, y)}{\Delta y_3} = \dfrac{\partial u}{\partial y} \\
& \dfrac{u(x, y + \Delta y_3) + u(x, y - \Delta y_3) - 2 u(x, y)}{\Delta y_3^2} = \dfrac{\partial^2 u}{\partial y^2}
\end{aligned}
\]

Therefore 1
\[
\begin{aligned}
\dfrac{u(x + \Delta x_3, y) - u(x, y)}{\Delta x_3} + \dfrac{u(x - \Delta x_2, y) - u(x, y)}{\Delta x_2} = - \beta\dfrac{\partial u}{\partial y} + \dfrac{\Delta x_2 + \Delta x_3}{2} \dfrac{\partial^2 u}{\partial x^2} + \dfrac{\Delta x_3}{2} \beta \left( 2\dfrac{\partial^2 u}{\partial x \partial y} + \beta\dfrac{\partial^2 u}{\partial y^2} \right)
\end{aligned}
\]

We could see that

\[
\begin{aligned}
\dfrac{\partial^2 u}{\partial x \partial y} + \beta\dfrac{\partial^2 u}{\partial y^2} = \dfrac{\partial}{\partial y}\left( \dfrac{\partial u}{\partial x} + \beta\dfrac{\partial u}{\partial y} \right) = & \dfrac{1}{\Delta y_3} \left[ \left(\dfrac{\partial u}{\partial x} + \beta\dfrac{\partial u}{\partial y} \right)(x, y) - \left(\dfrac{\partial u}{\partial x} + \beta\dfrac{\partial u}{\partial y} \right)(x, y - \Delta y_3) \right] \\
= & - \dfrac{u(x - \Delta x_2, y) - u(x, y)}{\Delta x_2 \Delta y_3} + \dfrac{u(x - \Delta x_2, y + \Delta y_3) - u(x, y + \Delta y_3)}{\Delta x_2 \Delta y_3}
\end{aligned}
\]

Therefore 2
\[
\begin{aligned}
\dfrac{u(x + \Delta x_3, y) - u(x, y)}{\Delta x_3} + \dfrac{u(x - \Delta x_2, y) - u(x, y)}{\Delta x_2} = & -\beta\dfrac{\partial u}{\partial y} + \dfrac{\Delta x_1 + \Delta x_2}{2} \dfrac{\partial^2 u}{\partial x^2} + \dfrac{\Delta x_2}{2} \beta \left( 2\dfrac{\partial^2 u}{\partial x \partial y} + 2 \beta\dfrac{\partial^2 u}{\partial y^2} \right) - \dfrac{\Delta x_3}{2} \beta^2 \dfrac{\partial^2 u}{\partial y^2} \\
= & - \beta \dfrac{u(x, y + \Delta y_3) - u(x, y)}{\Delta y_3} + \dfrac{\Delta x_2 + \Delta x_2}{2} \dfrac{\partial^2 u}{\partial x^2} \\
& + \Delta x_2 \beta \left( \dfrac{u(x - \Delta x_2, y + \Delta y_3) - u(x, y + \Delta y_3)}{\Delta x_2 \Delta y_3} - \dfrac{u(x - \Delta x_2, y) - u(x, y)}{\Delta x_2 \Delta y_3} \right) - \dfrac{\Delta x_2}{2} \beta^2 \dfrac{\partial^2 u}{\partial y^2} \\
\end{aligned}
\]

Therefore 3
\[
\begin{aligned}
2 \dfrac{u(x + \Delta x_3, y) - u(x, y)}{\Delta x_3 (\Delta x_2 + \Delta x_3)} + 2 \dfrac{u(x - \Delta x_2, y) - u(x, y)}{\Delta x_2 (\Delta x_2 + \Delta x_3)} - 2\beta\dfrac{u(x, y + \Delta y_3) - u(x, y)}{\Delta y_3 (\Delta x_1 + \Delta x_2)} - 2 \beta \left( \dfrac{u(x - \Delta x_2, y + \Delta y_3) - u(x, y + \Delta y_3)}{\Delta y_3 (\Delta x_1 + \Delta x_2)} - \dfrac{u(x - \Delta x_2, y) - u(x, y)}{\Delta y_3 (\Delta x_1 + \Delta x_2)} \right) = & \dfrac{\partial^2 u}{\partial x^2} - \dfrac{\Delta x_2}{\Delta x_2 + \Delta x_3} \beta^2 \dfrac{\partial^2 u}{\partial y^2}
\end{aligned}
\]

Therefore 4
\[
\begin{aligned}
& 2 \dfrac{u(x + \Delta x_3, y) - u(x, y)}{\Delta x_3 (\Delta x_2 + \Delta x_3)} + 2 \dfrac{u(x - \Delta x_2, y) - u(x, y)}{\Delta x_2 (\Delta x_2 + \Delta x_3)} - 2\beta\dfrac{u(x, y + \Delta y_3) - u(x, y)}{\Delta y_3 (\Delta x_2 + \Delta x_3)} - 2 \beta \left( \dfrac{u(x - \Delta x_2, y + \Delta y_3) - u(x, y + \Delta y_3)}{\Delta y_3 (\Delta x_2 + \Delta x_3)} - \dfrac{u(x - \Delta x_2, y) - u(x, y)}{\Delta y_3 (\Delta x_2 + \Delta x_3)} \right) \\
& \left( \dfrac{\Delta x_2}{\Delta x_2 + \Delta x_3}\beta^2 + 1 \right) \dfrac{u(x, y + \Delta y_3) + u(x, y - \Delta y_3) - 2 u(x, y)}{\Delta y_3^2}
= \dfrac{\partial^2 u}{\partial x^2} + \dfrac{\partial^2 u}{\partial y^2} = \Delta u
\end{aligned}
\]


===================================================================
===================================================================

We now express explicitly the discretization process in term of node value:

On the first rectangular domain: we define $\boxed{\beta_1 := \Delta x_1 / \Delta y_1}$ and $\boxed{\gamma_1 := -2(1 + \beta_1^2)}$. Then we have:
\[
\Delta \psi(x,y) = \dfrac{u(x + \Delta x_1, y) + u(x - \Delta x_1, y) - 2u(x, y)}{\Delta x_1^2} + \dfrac{u(x, y + \Delta y_1) + u(x, y - \Delta y_1) - 2u(x, y)}{\Delta y_1^2} = 0 \Longrightarrow -2(1+\beta_1^2)u(i,j) + \beta_1^2(u(i, j + 1) + u(i, j - 1)) + u(i+1,j) + u(i-1,j) = 0.
\]
\[
\boxed{
\gamma_1 u(i,j) + \beta_1^2(u(i, j + 1) + u(i, j - 1)) + u(i+1,j) + u(i-1,j) = 0
}.
\]

\[
\begin{aligned}
& \gamma_1 u(1, 1) + \beta_1^2(u(1, 2) + u(1, 0)) + u(2, 1) + u(0, 1) = 0\\
& \gamma_1 u(1, 2) + \beta_1^2(u(1, 3) + u(1, 1)) + u(2, 2) + u(0, 2) = 0\\
& \gamma_1 u(1, 3) + \beta_1^2(u(1, 4) + u(1, 2)) + u(2, 3) + u(0, 3) = 0.
\end{aligned}
\]
Therefore, for the rectangular domain:
\[
\begin{aligned}
& \gamma_1 u(1, 1) + u(2, 1) + \beta_1^2 u(1, 2) = - u(0, 1) - \beta_1^2 u(1, 0) \\
& \beta_1^2 u(1, 3) + \gamma_1 u(1, 2) + u(2, 2) + \beta_1^2 u(1, 1) = - u(0, 2) \\
& \beta_1^2 u(1, 2) + \gamma_1 u(1, 3) + u(2, 3) = - u(0, 3) - \beta_1^2 u(1, 4)
\end{aligned}.
\]

=====

On the oblique domain: we define $\boxed{\beta_2 := \Delta x_2 / \Delta y_2}$ and $\boxed{\gamma_2 := -2(1 + \beta \beta_2 + (\beta^2 + 1)\beta_2^2)}$, and $\boxed{\lambda := 2\beta\beta_2 + (\beta^2 + 1)\beta_2^2}$. Then we have:

\[
\begin{aligned}
& \dfrac{u(x + \Delta x_2, y) + u(x - \Delta x_2, y) - 2u(x,y)}{\Delta x_2^2} - 2\beta \dfrac{u(x + \Delta x_2, y + \Delta y_2) - u(x, y + \Delta y_2) - u(x + \Delta x_2, y) + u(x,y)}{\Delta x_2 \Delta y_2} + (\beta^2 + 1) \dfrac{u(x, y + \Delta y_2) + u(x, y - \Delta y_2) - 2u(x,y)}{\Delta y_2^2} = 0 \\
& \dfrac{u(x + \Delta x_2, y) + u(x - \Delta x_2, y) - 2u(x,y)}{\beta_2^2 \Delta y_2^2} - 2\beta \dfrac{u(x + \Delta x_2, y + \Delta y_2) - u(x, y + \Delta y_2) - u(x + \Delta x_2, y) + u(x,y)}{\beta_2 \Delta y_2 \Delta y_2} + (\beta^2 + 1) \dfrac{u(x, y + \Delta y_2) + u(x, y - \Delta y_2) - 2u(x,y)}{\Delta y_2^2} = 0 \\
& u(x + \Delta x_2, y) + u(x - \Delta x_2, y) - 2u(x,y) - 2\beta\beta_2 \left(u(x + \Delta x_2, y + \Delta y_2) - u(x, y + \Delta y_2) - u(x + \Delta x_2, y) + u(x,y) \right) + (\beta^2 + 1)\beta_2^2 \left( u(x, y + \Delta y_2) + u(x, y - \Delta y_2) - 2u(x,y) \right) = 0 \\
\end{aligned}
\]
In term of nodal values, we have
\[
\begin{aligned}
& u(i + 1, j) + u(i - 1, j) - 2u(i, j) - 2\beta\beta_2 \left(u(i + 1, j + 1) - u(i, j + 1) - u(i + 1, j) + u(i, j) \right) + (\beta^2 + 1)\beta_2^2 \left( u(i, j + 1) + u(i, j - 1) - 2u(i, j) \right) = 0 \\
& -2(1 + \beta \beta_2 + (\beta^2 + 1)\beta_2^2) u(i, j) + (1 + 2\beta\beta_2) u(i + 1, j) + u(i - 1, j) + (2\beta\beta_2 + (\beta^2 + 1)\beta_2^2) u(i, j + 1) - 2\beta\beta_2 u(i + 1, j + 1) + (\beta^2 + 1)\beta_2^2 u(i, j - 1) = 0 \\
& u(i - 1, j) -2(1 + \beta \beta_2 + (\beta^2 + 1)\beta_2^2) u(i, j) + (1 + 2\beta\beta_2) u(i + 1, j) + (2\beta\beta_2 + (\beta^2 + 1)\beta_2^2) u(i, j + 1) - 2\beta\beta_2 u(i + 1, j + 1) + (\beta^2 + 1)\beta_2^2 u(i, j - 1) = 0 \\
\end{aligned}
\]
Therefore
\[
\begin{aligned}
u(i - 1, j) + \gamma_2 u(i, j) + (1 + 2\beta\beta_2) u(i + 1, j) + \lambda u(i, j + 1) - 2\beta\beta_2 u(i + 1, j + 1) + (\beta^2 + 1)\beta_2^2 u(i, j - 1) = 0
\end{aligned}
\]

Specifically,
\[
\begin{aligned}
& u(2, 1) + \gamma_2 u(3, 1) + (1 + 2\beta\beta_2) u(4, 1) + \lambda u(3, 2) - 2\beta\beta_2 u(4, 2) + (\beta^2 + 1)\beta_2^2 u(3, 0) = 0 \\
& u(2, 2) + \gamma_2 u(3, 2) + (1 + 2\beta\beta_2) u(4, 2) + \lambda u(3, 3) - 2\beta\beta_2 u(4, 3) + (\beta^2 + 1)\beta_2^2 u(3, 1) = 0 \\
& u(2, 3) + \gamma_2 u(3, 3) + (1 + 2\beta\beta_2) u(4, 3) + \lambda u(3, 4) - 2\beta\beta_2 u(4, 4) + (\beta^2 + 1)\beta_2^2 u(3, 2) = 0 \\
& u(3, 1) + \gamma_2 u(4, 1) + (1 + 2\beta\beta_2) u(5, 1) + \lambda u(4, 2) - 2\beta\beta_2 u(5, 2) + (\beta^2 + 1)\beta_2^2 u(4, 0) = 0 \\
& u(3, 2) + \gamma_2 u(4, 2) + (1 + 2\beta\beta_2) u(5, 2) + \lambda u(4, 3) - 2\beta\beta_2 u(5, 3) + (\beta^2 + 1)\beta_2^2 u(4, 1) = 0 \\
& u(3, 3) + \gamma_2 u(4, 3) + (1 + 2\beta\beta_2) u(5, 3) + \lambda u(4, 4) - 2\beta\beta_2 u(5, 4) + (\beta^2 + 1)\beta_2^2 u(4, 2) = 0 \\
\end{aligned}
\]

\[
\begin{aligned}
& u(2, 1) + \gamma_2 u(3, 1) + (1 + 2\beta\beta_2) u(4, 1) + \lambda u(3, 2) - 2\beta\beta_2 u(4, 2) = -(\beta^2 + 1)\beta_2^2 u(3, 0) \\
& u(2, 2) + \gamma_2 u(3, 2) + (1 + 2\beta\beta_2) u(4, 2) + \lambda u(3, 3) - 2\beta\beta_2 u(4, 3) + (\beta^2 + 1)\beta_2^2 u(3, 1) = 0 \\
& u(2, 3) + \gamma_2 u(3, 3) + (1 + 2\beta\beta_2) u(4, 3) + (\beta^2 + 1)\beta_2^2 u(3, 2) = - \lambda u(3, 4) + 2\beta\beta_2 u(4, 4) \\
& u(3, 1) + \gamma_2 u(4, 1) + (1 + 2\beta\beta_2) u(5, 1) + \lambda u(4, 2) - 2\beta\beta_2 u(5, 2) = -(\beta^2 + 1)\beta_2^2 u(4, 0) \\
& u(3, 2) + \gamma_2 u(4, 2) + (1 + 2\beta\beta_2) u(5, 2) + \lambda u(4, 3) - 2\beta\beta_2 u(5, 3) + (\beta^2 + 1)\beta_2^2 u(4, 1) = 0 \\
& u(3, 3) + \gamma_2 u(4, 3) + (1 + 2\beta\beta_2) u(5, 3) + (\beta^2 + 1)\beta_2^2 u(4, 2) = - \lambda u(4, 4) + 2\beta\beta_2 u(5, 4) \\
\end{aligned}
\]

===
On the second rectangular domain: we define $\boxed{\beta_3 := \Delta x_3 / \Delta y_3}$ and $\boxed{\gamma_3 := -2(1 + \beta_3^2)}$. Then we have:
\[
\boxed{
\gamma_3 u(i,j) + \beta_3^2(u(i, j + 1) + u(i, j - 1)) + u(i+1,j) + u(i-1,j) = 0
}.
\]

\[
\begin{aligned}
& \gamma_3 u(6, 1) + \beta_3^2(u(6, 2) + u(6, 0)) + u(7, 1) + u(5, 1) = 0\\
& \gamma_3 u(6, 2) + \beta_3^2(u(6, 3) + u(6, 1)) + u(7, 2) + u(5, 2) = 0\\
& \gamma_3 u(6, 3) + \beta_3^2(u(6, 4) + u(6, 2)) + u(7, 3) + u(5, 3) = 0.
\end{aligned}
\]
Therefore, for the rectangular domain:
\[
\begin{aligned}
& \gamma_3 u(6, 1) + \beta_3^2 u(6, 2) + u(5, 1) = - \beta_3^2 u(6, 0) - u(7, 1)\\
& \gamma_3 u(6, 2) + \beta_3^2 u(6, 3) + \beta_3^2 u(6, 1) + u(5, 2) = - u(7, 2)\\
& \gamma_3 u(6, 3) + \beta_3^2 u(6, 2) + u(5, 3) = -\beta_3^2 u(6, 4) - u(7, 3).
\end{aligned}
\]

==========

On the first interface $\boxed{k_1 := \Delta x_2/ \Delta x_1}$:
\[
\begin{aligned}
& 2 \dfrac{u(x + \Delta x_2, y) - u(x, y)}{k_1 \beta_1 \Delta y_1 (1 + k_1) \beta_1 \Delta y_1} + 2 \dfrac{u(x - \Delta x_1, y) - u(x, y)}{\beta_1\Delta y_1 (1 + k_1)\beta_1\Delta y_1} - 2\beta\dfrac{u(x, y + \Delta y_1) - u(x, y)}{\Delta y_1 (1 + k_1) \beta_1\Delta y_1} - 2 \beta \dfrac{u(x + \Delta x_2, y + \Delta y_1) - u(x, y + \Delta y_1) - u(x + \Delta x_2, y) + u(x, y)}{\Delta y_1 (1 + k_1)\beta_1 \Delta y_1} \\
& + \left( \dfrac{k_1}{1 + k_1} \beta^2 + 1 \right) \dfrac{u(x, y + \Delta y_1) + u(x, y - \Delta y_1) - 2 u(x, y)}{\Delta y_1^2} = 0
\end{aligned}
\]

Therefore
\[
\begin{aligned}
& 2 \dfrac{u(x + \Delta x_2, y) - u(x, y)}{k_1 (1 + k_1) \beta_1^2} + 2 \dfrac{u(x - \Delta x_1, y) - u(x, y)}{(1 + k_1) \beta_1^2} - 2\beta\dfrac{u(x, y + \Delta y_1) - u(x, y)}{(1 + k_1) \beta_1} - 2 \beta \dfrac{u(x + \Delta x_2, y + \Delta y_1) - u(x, y + \Delta y_1) - u(x + \Delta x_2, y) + u(x, y)}{(1 + k_1)\beta_1} \\
& + \left( \dfrac{k_1}{1 + k_1} \beta^2 + 1 \right) \left(u(x, y + \Delta y_1) + u(x, y - \Delta y_1) - 2 u(x, y) \right) = 0
\end{aligned}
\]

\[
\begin{aligned}
& -2[1 + k_1 + k_1 \beta_1^2 (1 + k_1 + k_1 \beta^2)] u(x,y) + 2(1 + k_1 \beta_1 \beta) u(x + \Delta x_2, y) + 2k_1 u(x - \Delta x_1, y) + k_1 \beta_1^2 (1 + k_1 + k_1\beta^2) u(x, y + \Delta y_1) + k_1 \beta_1^2 (1 + k_1 + k_1\beta^2) u(x, y - \Delta y_1) - 2 k_1 \beta_1 \beta u(x + \Delta x_2, y + \Delta y_1) = 0.
\end{aligned}
\]

By defining $\boxed{\gamma_{12} := -2[1 + k_1 + k_1 \beta_1^2 (1 + k_1 + k_1 \beta^2)]}$, $\boxed{\beta_{12} := k_1 \beta_1^2 (1 + k_1 + k_1\beta^2)}$

\[
\begin{aligned}
& \gamma_{12} u(x,y) + 2(1 + k_1 \beta_1 \beta) u(x + \Delta x_2, y) + 2k_1 u(x - \Delta x_1, y) + \beta_{12} u(x, y + \Delta y_1) + \beta_{12} u(x, y - \Delta y_1) - 2 k_1 \beta_1 \beta u(x + \Delta x_2, y + \Delta y_1) = 0.
\end{aligned}
\]

This gives us:
\[
\begin{aligned}
& \gamma_{12} u(i, j) + 2(1 + k_1 \beta_1 \beta) u(i + 1, j) + 2k_1 u(i - 1, j) + \beta_{12} u(i, j + 1) + \beta_{12} u(i, j - 1) - 2 k_1 \beta_1 \beta u(i + 1, j + 1) = 0.
\end{aligned}
\]

Specifically
\[
\begin{aligned}
& \gamma_{12} u(2, 1) + 2(1 + k_1 \beta_1 \beta) u(3, 1) + 2k_1 u(1, 1) + \beta_{12} u(2, 2) + \beta_{12} u(2, 0) - 2 k_1 \beta_1 \beta u(3, 2) = 0\\
& \gamma_{12} u(2, 2) + 2(1 + k_1 \beta_1 \beta) u(3, 2) + 2k_1 u(1, 2) + \beta_{12} u(2, 3) + \beta_{12} u(2, 1) - 2 k_1 \beta_1 \beta u(3, 3) = 0\\
& \gamma_{12} u(2, 3) + 2(1 + k_1 \beta_1 \beta) u(3, 3) + 2k_1 u(1, 3) + \beta_{12} u(2, 4) + \beta_{12} u(2, 2) - 2 k_1 \beta_1 \beta u(3, 4) = 0.
\end{aligned}
\]

\[
\begin{aligned}
& \gamma_{12} u(2, 1) + 2(1 + k_1 \beta_1 \beta) u(3, 1) + 2k_1 u(1, 1) + \beta_{12} u(2, 2) - 2 k_1 \beta_1 \beta u(3, 2) = - \beta_{12} u(2, 0)\\
& \gamma_{12} u(2, 2) + 2(1 + k_1 \beta_1 \beta) u(3, 2) + 2k_1 u(1, 2) + \beta_{12} u(2, 3) + \beta_{12} u(2, 1) - 2 k_1 \beta_1 \beta u(3, 3) = 0\\
& \gamma_{12} u(2, 3) + 2(1 + k_1 \beta_1 \beta) u(3, 3) + 2k_1 u(1, 3) + \beta_{12} u(2, 2) = - \beta_{12} u(2, 4) + 2 k_1 \beta_1 \beta u(3, 4).
\end{aligned}
\]


==========

On the second interface $\boxed{k_2 := \Delta x_3/ \Delta x_2}$:

\[
\begin{aligned}
& 2 \dfrac{u(x + \Delta x_3, y) - u(x, y)}{\Delta y_3 \beta_3 (\frac{1}{k_2} + 1) \Delta y_3 \beta_3} + 2 \dfrac{u(x - \Delta x_2, y) - u(x, y)}{\frac{\beta_3}{k_2}\Delta y_3 (\frac{1}{k_2} + 1) \Delta y_3 \beta_3} - 2\beta\dfrac{u(x, y + \Delta y_3) - u(x, y)}{\Delta y_3 (\frac{1}{k_2} + 1) \Delta y_3 \beta_3} - 2 \beta \dfrac{u(x - \Delta x_2, y + \Delta y_3) - u(x, y + \Delta y_3) - u(x - \Delta x_2, y) + u(x, y)}{\Delta y_3 (\frac{1}{k_2} + 1) \Delta y_3 \beta_3} \\
& \left( \dfrac{\beta^2}{1 + k_2} + 1 \right) \dfrac{u(x, y + \Delta y_3) + u(x, y - \Delta y_3) - 2 u(x, y)}{\Delta y_3^2}
= \dfrac{\partial^2 u}{\partial x^2} + \dfrac{\partial^2 u}{\partial y^2} = \Delta u
\end{aligned}
\]


\[
\begin{aligned}
& 2 \dfrac{u(x + \Delta x_3, y) - u(x, y)}{\beta_3 (\frac{1}{k_2} + 1)\beta_3} + 2 \dfrac{u(x - \Delta x_2, y) - u(x, y)}{\frac{\beta_3}{k_2} (\frac{1}{k_2} + 1)\beta_3} - 2\beta\dfrac{u(x, y + \Delta y_3) - u(x, y)}{(\frac{1}{k_2} + 1) \beta_3} - 2 \beta \dfrac{u(x - \Delta x_2, y + \Delta y_3) - u(x, y + \Delta y_3) - u(x - \Delta x_2, y) + u(x, y)}{(\frac{1}{k_2} + 1) \beta_3} \\
& + \left( \dfrac{\beta^2}{1 + k_2} + 1 \right) \left(u(x, y + \Delta y_3) + u(x, y - \Delta y_3) - 2 u(x, y) \right) = 0
\end{aligned}
\]

\[
\begin{aligned}
& 2 \dfrac{u(x + \Delta x_3, y) - u(x, y)}{\beta_3 (\frac{1}{k_2} + 1)\beta_3} + 2 \dfrac{u(x - \Delta x_2, y) - u(x, y)}{\frac{\beta_3}{k_2} (\frac{1}{k_2} + 1)\beta_3} - 2 \beta \dfrac{u(x - \Delta x_2, y + \Delta y_3) - u(x - \Delta x_2, y)}{(\frac{1}{k_2} + 1) \beta_3} \\
& + \left( \dfrac{\beta^2}{1 + k_2} + 1 \right) \left(u(x, y + \Delta y_3) + u(x, y - \Delta y_3) - 2 u(x, y) \right) = 0
\end{aligned}
\]

\[
\begin{aligned}
-2[k_2 + k_2^2 + \beta_3^2(1 + k_2 + \beta^2)] u(x, y) + 2k_2 u(x + \Delta x_3, y) + 2k_2(k_2 + \beta\beta_3) u(x - \Delta x_2, y) + (1 + k_2 + \beta^2) \beta_3^2 u(x, y + \Delta y_3) + (1 + k_2 + \beta^2) \beta_3^2 u(x, y - \Delta y_3) + 2k_2\beta_3\beta u(x - \Delta x_2, y + \Delta y_3)= 0.
\end{aligned}
\]
Be defining $\boxed{\gamma_{23} := -2[k_2 + k_2^2 + \beta_3^2(1 + k_2 + \beta^2)]}$, $\boxed{\beta_{23} := (1 + k_2 + \beta^2) \beta_3^2}$

\[
\begin{aligned}
\gamma_{23} u(x, y) + 2k_2 u(x + \Delta x_3, y) + 2k_2(k_2 + \beta\beta_3) u(x - \Delta x_2, y) + \beta_{23} u(x, y + \Delta y_3) + \beta_{23} \beta_3^2 u(x, y - \Delta y_3) + 2k_2\beta_3\beta u(x - \Delta x_2, y + \Delta y_3)= 0.
\end{aligned}
\]

This gives us
\[
\begin{aligned}
\gamma_{23} u(i, j) + 2k_2 u(i + 1, j) + 2k_2(k_2 + \beta\beta_3) u(i - 1, j) + \beta_{23} u(i, j + 1) + \beta_{23} \beta_3^2 u(i, j - 1) + 2k_2\beta_3\beta u(i - 1, j + 1) = 0.
\end{aligned}
\]

Specifically

\[
\begin{aligned}
& \gamma_{23} u(5, 1) + 2k_2 u(6, 1) + 2k_2(k_2 + \beta\beta_3) u(4, 1) + \beta_{23} u(5, 2) + \beta_{23} \beta_3^2 u(5, 0) + 2k_2\beta_3\beta u(4, 2) = 0 \\
& \gamma_{23} u(5, 2) + 2k_2 u(6, 2) + 2k_2(k_2 + \beta\beta_3) u(4, 2) + \beta_{23} u(5, 3) + \beta_{23} \beta_3^2 u(5, 1) + 2k_2\beta_3\beta u(4, 3) = 0 \\
& \gamma_{23} u(5, 3) + 2k_2 u(6, 3) + 2k_2(k_2 + \beta\beta_3) u(4, 3) + \beta_{23} u(5, 4) + \beta_{23} \beta_3^2 u(5, 2) + 2k_2\beta_3\beta u(4, 4) = 0 \\
\end{aligned}
\]

\[
\begin{aligned}
& \gamma_{23} u(5, 1) + 2k_2 u(6, 1) + 2k_2(k_2 + \beta\beta_3) u(4, 1) + \beta_{23} u(5, 2) + 2k_2\beta_3\beta u(4, 2) = - \beta_{23} \beta_3^2 u(5, 0) \\
& \gamma_{23} u(5, 2) + 2k_2 u(6, 2) + 2k_2(k_2 + \beta\beta_3) u(4, 2) + \beta_{23} u(5, 3) + \beta_{23} \beta_3^2 u(5, 1) + 2k_2\beta_3\beta u(4, 3) = 0 \\
& \gamma_{23} u(5, 3) + 2k_2 u(6, 3) + 2k_2(k_2 + \beta\beta_3) u(4, 3) + \beta_{23} \beta_3^2 u(5, 2) = - \beta_{23} u(5, 4) - 2k_2\beta_3\beta u(4, 4) = 0 \\
\end{aligned}
\]

===================================================================
===================================================================

For the oblique grid exclusively:

\[
\begin{bmatrix}
\gamma_1 & 1 & 0 & 0 & 0 & 0 & \beta_1^2 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
1 & \gamma_1 & 1 & 0 & 0 & 0 & 0 & \beta_1^2 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 1 & \gamma_2 & 1 + \beta\beta_2 & 0 & 0 & 0 & 0 & \lambda & -2\beta\beta_2 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 1 & \gamma_1 & 1 & 0 & 0 & 0 & 0 & \beta_1^2 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 1 & \gamma_1 & 1 & 0 & 0 & 0 & 0 & \beta_1^2 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & \gamma_3 & 1 & 0 & 0 & 0 & 0 & \beta_3^2 & 0 & 0 & 0 & 0 & 0 & 0 \\
\beta_1^2 & 0 & 0 & 0 & 0 & 1 & \gamma_1 & 1 & 0 & 0 & 0 & 0 & \beta_1^2 & 0 & 0 & 0 & 0 & 0 \\
0 & \beta_3^2 & 0 & 0 & 0 & 0 & 1 & \gamma_1 & 1 & 0 & 0 & 0 & 0 & \beta_1^2 & 0 & 0 & 0 & 0 \\
0 & 0 & (\beta^2 + 1)\beta_2^2 & 0 & 0 & 0 & 0 & 1 & \gamma_2 & 1 + 2\beta\beta_2 & 0 & 0 & 0 & 0 & \lambda & -2\beta\beta_2 & 0 & 0 \\
0 & 0 & 0 & \beta_3^2 & 0 & 0 & 0 & 0 & 1 & \gamma_2 & 1 & 0 & 0 & 0 & 0 & \beta_1^2 & 0 & 0 \\
0 & 0 & 0 & 0 & \beta_3^2 & 0 & 0 & 0 & 0 & 1 & \gamma_1 & 1 & 0 & 0 & 0 & 0 & \beta_1^2 & 0 \\
0 & 0 & 0 & 0 & 0 & \beta_3^2 & 0 & 0 & 0 & 0 & 1 & \gamma_3 & 0 & 0 & 0 & 0 & 0 & \beta_3^2 \\
0 & 0 & 0 & 0 & 0 & 0 & \beta_1^2 & 0 & 0 & 0 & 0 & 1 & \gamma_1 & 1 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & \beta_3^2 & 0 & 0 & 0 & 0 & 1 & \gamma_1 & 1 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & (\beta^2 + 1)\beta_2^2 & 0 & 0 & 0 & 0 & 1 & \gamma_2 & 1 + 2\beta\beta_2 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & (\beta^2 + 1)\beta_2^2 & 0 & 0 & 0 & 0 & 1 & \gamma_2 & 1 + 2\beta\beta_2 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \beta_3^2 & 0 & 0 & 0 & 0 & 1 & \gamma_1 & 1 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \beta_3^2 & 0 & 0 & 0 & 0 & 1 & \gamma_3
\end{bmatrix}\begin{bmatrix}
u_{1,1} \\
u_{2,1} \\
u_{3,1} \\
u_{4,1} \\
u_{5,1} \\
u_{6,1} \\
u_{1,2} \\
u_{2,2} \\
u_{3,2} \\
u_{4,2} \\
u_{5,2} \\
u_{6,2} \\
u_{1,3} \\
u_{2,3} \\
u_{3,3} \\
u_{4,3} \\
u_{5,3} \\
u_{6,3} \\
\end{bmatrix}=\begin{bmatrix}
-u_{0,1} - \beta_1^2 u_{1,0} \\
u_{2,1} \\
-(\beta^2 + 1)\beta_2^2 u_{3,0} \\
u_{4,1} \\
u_{5,1} \\
-u_{7,1} - \beta_3^2 u_{6,0}\\
-u_{0,2} \\
u_{2,2} \\
0 \\
u_{4,2} \\
u_{5,2} \\
-u_{7,2} \\
-u_{0,3} - \beta_1^2 u_{1,4}\\
u_{2,3} \\
-2\beta\beta_2 u_{4,4} - \lambda u_{3,4} \\
-\lambda u_{4,4} + 2\beta\beta_2 u_{5,4} \\
u_{5,3} \\
-u_{7,3} - \beta_3^2 u_{6,4} \\
\end{bmatrix}
\]




\[
\begin{bmatrix}
\gamma_1 & 1 & 0 & 0 & 0 & 0 & \beta_1^2 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
2k_1 & \gamma_{12} & 2(1 + k_1\beta_1\beta) & 0 & 0 & 0 & 0 & \beta_{12} & -2k_1\beta_1\beta & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 1 & \gamma_2 & (1 + 2\beta\beta_2) & 0 & 0 & 0 & 0 & \lambda & -2\beta\beta_2 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 1 & \gamma_2 & 1 + 2\beta\beta_2 & 0 & 0 & 0 & 0 & \lambda & -2\beta\beta_2 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 2k_2(k_2 + \beta\beta_3) & \gamma_{23} & 2k_2 & 0 & 0 & 0 & 2k_2\beta_3\beta & \beta_{23} & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & \gamma_3 & 0 & 0 & 0 & 0 & 0 & \beta_3^2 & 0 & 0 & 0 & 0 & 0 & 0 \\
\beta_1^2 & 0 & 0 & 0 & 0 & 0 & \gamma_1 & 1 & 0 & 0 & 0 & 0 & 0 & \beta_1^2 & 0 & 0 & 0 & 0 \\
0 & \beta_{12} & 0 & 0 & 0 & 0 & 2k_1 & \gamma_{12} & 2(1 + k_1\beta_1\beta) & 0 & 0 & 0 & 0 & \beta_{12} & -2k_1\beta_1\beta & 0 & 0 & 0 \\
0 & 0 & (\beta^2 + 1)\beta_2^2 & 0 & 0 & 0 & 0 & 1 & \gamma_2 & (1 + 2\beta\beta_2) & 0 & 0 & 0 & 0 & \lambda & -2\beta\beta_2 & 0 & 0 \\
0 & 0 & 0 & (\beta^2 + 1)\beta_2^2 & 0 & 0 & 0 & 0 & 1 & \gamma_2 & (1 + 2\beta\beta_2) & 0 & 0 & 0 & 0 & \lambda & -2\beta\beta_2 & 0 \\
0 & 0 & 0 & 0 & \beta_{23}\beta_3^2 & 0 & 0 & 0 & 0 & 2k_2(k_2 + \beta\beta_3) & \gamma_{23} & 2k_2 & 0 & 0 & 0 & 2k_2\beta_3\beta & \beta_{23} & 0 \\
0 & 0 & 0 & 0 & 0 & \beta_3^2 & 0 & 0 & 0 & 0 & 1 & \gamma_3 & 0 & 0 & 0 & 0 & 0 & \beta_3^2 \\
0 & 0 & 0 & 0 & 0 & 0 & \beta_1^2 & 0 & 0 & 0 & 0 & 0 & \gamma_1 & 1 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \beta_{12} & 0 & 0 & 0 & 0 & 2k_1 & \gamma_{12} & 2(1 + k_1\beta_1\beta) & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & (\beta^2 + 1)\beta_2^2 & 0 & 0 & 0 & 0 & 1 & \gamma_2 & (1 + 2\beta\beta_2) & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & (\beta^2 + 1)\beta_2^2 & 0 & 0 & 0 & 0 & 1 & \gamma_2 & (1 + 2\beta\beta_2) & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \beta_{23}\beta_3^2 & 0 & 0 & 0 & 0 & 2k_2(k_2 + \beta\beta_3) & \gamma_{23} & 2k_2 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \beta_3^2 & 0 & 0 & 0 & 0 & 1 & \gamma_3
\end{bmatrix}\begin{bmatrix}
u_{1,1} \\
u_{2,1} \\
u_{3,1} \\
u_{4,1} \\
u_{5,1} \\
u_{6,1} \\
u_{1,2} \\
u_{2,2} \\
u_{3,2} \\
u_{4,2} \\
u_{5,2} \\
u_{6,2} \\
u_{1,3} \\
u_{2,3} \\
u_{3,3} \\
u_{4,3} \\
u_{5,3} \\
u_{6,3}
\end{bmatrix}=\begin{bmatrix}
-u_{0,1} - \beta_1^2 u_{1,0} \\
-\beta_{12} u_{2,0} \\
-(\beta^2 + 1)\beta_2^2 u_{3,0} \\
-(\beta^2 + 1)\beta_2^2 u_{4,0} \\
-\beta_{23}\beta_3^2 u_{5,0} \\
-\beta_3^2 u_{6,0} - u_{7,1} \\
u_{0,2} \\
0 \\
0 \\
0 \\
0 \\
u_{7,2} \\
-u_{0,3} - \beta_1^2 u_{1,4} \\
-\beta_{12} u_{2,4} + 2k_1\beta_1\beta u_{3,4} \\
-\lambda u_{3,4} + 2\beta\beta_2 u_{4,4} \\
-\lambda u_{4,4} + 2\beta\beta_2 u_{5,4} \\
-\beta_{23} u_{5,4} - 2k_2\beta_3\beta u_{4,4} \\
-\beta_3^2 u_{6,4} - u_{7,3}
\end{bmatrix}.
\]



