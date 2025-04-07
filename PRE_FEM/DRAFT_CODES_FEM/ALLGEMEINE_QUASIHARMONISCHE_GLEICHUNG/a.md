\[
\begin{bmatrix}
    \partial_x N_1 & \partial_y N_1 \\
    \partial_x N_2 & \partial_y N_2 \\
    \partial_x N_3 & \partial_y N_3 \\
    \partial_x N_4 & \partial_y N_4 \\
\end{bmatrix}
\begin{bmatrix}
    \partial_x N_1 & \partial_x N_2 & \partial_x N_3 & \partial_x N_4 \\
    \partial_y N_1 & \partial_y N_2 & \partial_y N_3 & \partial_y N_4 \\
\end{bmatrix}=\begin{bmatrix}
    \partial_x N_1 \partial_x N_1 + \partial_y N_1  \partial_y N_1 & \partial_x N_1 \partial_x N_2 + \partial_y N_1  \partial_y N_2 & \partial_x N_1 \partial_x N_3 + \partial_y N_1  \partial_y N_3 & \partial_x N_1 \partial_x N_4 + \partial_y N_1  \partial_y N_4 \\
    \partial_x N_2 \partial_x N_1 + \partial_y N_2  \partial_y N_1 & \partial_x N_2 \partial_x N_2 + \partial_y N_2  \partial_y N_2 & \partial_x N_2 \partial_x N_3 + \partial_y N_2  \partial_y N_3 & \partial_x N_2 \partial_x N_4 + \partial_y N_2  \partial_y N_4 \\
    \partial_x N_3 \partial_x N_1 + \partial_y N_3  \partial_y N_1 & \partial_x N_3 \partial_x N_2 + \partial_y N_3  \partial_y N_2 & \partial_x N_3 \partial_x N_3 + \partial_y N_3  \partial_y N_3 & \partial_x N_3 \partial_x N_4 + \partial_y N_3  \partial_y N_4 \\
    \partial_x N_4 \partial_x N_1 + \partial_y N_4  \partial_y N_1 & \partial_x N_4 \partial_x N_2 + \partial_y N_4  \partial_y N_2 & \partial_x N_4 \partial_x N_3 + \partial_y N_4  \partial_y N_3 & \partial_x N_4 \partial_x N_4 + \partial_y N_4  \partial_y N_4 \\
\end{bmatrix}
\]

For the local stiffness matrix, the integrand is
\[
\begin{align*}
& k_{ij} = (\partial_x N_i \partial_x N_j + \partial_y N_i \partial_y N_j)\left(\sum_{k = 1}^4 N_k(\xi, \eta)y_i\right) \\
& k_{11} = (\partial_x N_1 \partial_x N_1 + \partial_y N_1 \partial_y N_1)\left(N_1 y_1 + N_2 y_2 + N_3 y_3 + N_4 y_4\right) \\
& k_{12} = (\partial_x N_1 \partial_x N_2 + \partial_y N_1 \partial_y N_2)\left(N_1 y_1 + N_2 y_2 + N_3 y_3 + N_4 y_4\right) \\
& k_{13} = (\partial_x N_1 \partial_x N_3 + \partial_y N_1 \partial_y N_3)\left(N_1 y_1 + N_2 y_2 + N_3 y_3 + N_4 y_4\right) \\
& k_{14} = (\partial_x N_1 \partial_x N_4 + \partial_y N_1 \partial_y N_4)\left(N_1 y_1 + N_2 y_2 + N_3 y_3 + N_4 y_4\right) \\
& k_{21} = (\partial_x N_2 \partial_x N_1 + \partial_y N_2 \partial_y N_1)\left(N_1 y_1 + N_2 y_2 + N_3 y_3 + N_4 y_4\right) \\
& k_{22} = (\partial_x N_2 \partial_x N_2 + \partial_y N_2 \partial_y N_2)\left(N_1 y_1 + N_2 y_2 + N_3 y_3 + N_4 y_4\right) \\
& k_{23} = (\partial_x N_2 \partial_x N_3 + \partial_y N_2 \partial_y N_3)\left(N_1 y_1 + N_2 y_2 + N_3 y_3 + N_4 y_4\right) \\
& k_{24} = (\partial_x N_2 \partial_x N_4 + \partial_y N_2 \partial_y N_4)\left(N_1 y_1 + N_2 y_2 + N_3 y_3 + N_4 y_4\right) \\
& k_{31} = (\partial_x N_3 \partial_x N_1 + \partial_y N_3 \partial_y N_1)\left(N_1 y_1 + N_2 y_2 + N_3 y_3 + N_4 y_4\right) \\
& k_{32} = (\partial_x N_3 \partial_x N_2 + \partial_y N_3 \partial_y N_2)\left(N_1 y_1 + N_2 y_2 + N_3 y_3 + N_4 y_4\right) \\
& k_{33} = (\partial_x N_3 \partial_x N_3 + \partial_y N_3 \partial_y N_3)\left(N_1 y_1 + N_2 y_2 + N_3 y_3 + N_4 y_4\right) \\
& k_{34} = (\partial_x N_3 \partial_x N_4 + \partial_y N_3 \partial_y N_4)\left(N_1 y_1 + N_2 y_2 + N_3 y_3 + N_4 y_4\right) \\
& k_{41} = (\partial_x N_4 \partial_x N_1 + \partial_y N_4 \partial_y N_1)\left(N_1 y_1 + N_2 y_2 + N_3 y_3 + N_4 y_4\right) \\
& k_{42} = (\partial_x N_4 \partial_x N_2 + \partial_y N_4 \partial_y N_2)\left(N_1 y_1 + N_2 y_2 + N_3 y_3 + N_4 y_4\right) \\
& k_{43} = (\partial_x N_4 \partial_x N_3 + \partial_y N_4 \partial_y N_3)\left(N_1 y_1 + N_2 y_2 + N_3 y_3 + N_4 y_4\right) \\
& k_{44} = (\partial_x N_4 \partial_x N_4 + \partial_y N_4 \partial_y N_4)\left(N_1 y_1 + N_2 y_2 + N_3 y_3 + N_4 y_4\right) \\
\end{align*}
\]

\[
[k_{ij}] = \left(N_1 y_1 + N_2 y_2 + N_3 y_3 + N_4 y_4\right)\begin{bmatrix}
    \partial_x N_1 & \partial_y N_1 \\
    \partial_x N_2 & \partial_y N_2 \\
    \partial_x N_3 & \partial_y N_3 \\
    \partial_x N_4 & \partial_y N_4 \\
\end{bmatrix}
\begin{bmatrix}
    \partial_x N_1 & \partial_x N_2 & \partial_x N_3 & \partial_x N_4 \\
    \partial_y N_1 & \partial_y N_2 & \partial_y N_3 & \partial_y N_4 \\
\end{bmatrix} = AB^TB
\]
where
\[
A = N_1 y_1 + N_2 y_2 + N_3 y_3 + N_4 y_4 =
\begin{bmatrix}
N_1 & N_2 & N_3 & N_4
\end{bmatrix}
\begin{bmatrix}
y_1\\
y_2\\
y_3\\
y_4\\
\end{bmatrix}
\]
\[
B = \begin{bmatrix}
    \partial_x N_1 & \partial_x N_2 & \partial_x N_3 & \partial_x N_4 \\
    \partial_y N_1 & \partial_y N_2 & \partial_y N_3 & \partial_y N_4 \\
\end{bmatrix}=\begin{bmatrix}
    \partial_x \xi & \partial_x \eta \\
    \partial_y \xi & \partial_y \eta \\
\end{bmatrix}\begin{bmatrix}
    \partial_\xi N_1 & \partial_\xi N_2 & \partial_\xi N_3 & \partial_\xi N_4 \\
    \partial_\eta N_1 & \partial_\eta N_2 & \partial_\eta N_3 & \partial_\eta N_4 \\
\end{bmatrix}
\]

We have Jacobian matrix
\[
J = \begin{bmatrix}
    \partial_\xi x & \partial_\xi y \\
    \partial_\eta x & \partial_\eta y \\
\end{bmatrix} \quad
J^{-1} = \begin{bmatrix}
    \partial_x \xi & \partial_x \eta \\
    \partial_y \xi & \partial_y \eta \\
\end{bmatrix}
\]
and the Parent Gradient Matrix
\[
G = \begin{bmatrix}
    \partial_\xi N_1 & \partial_\xi N_2 & \partial_\xi N_3 & \partial_\xi N_4 \\
    \partial_\eta N_1 & \partial_\eta N_2 & \partial_\eta N_3 & \partial_\eta N_4 \\
\end{bmatrix}
\]

Therefore
\[
B = J^{-1}G
\]