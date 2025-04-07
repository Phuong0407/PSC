\[
\begin{aligned}
& \partial_t u_1 = -R \partial_1 (u_1^2) - R \partial_2 (u_1 u_2) - R \partial_2 (u_1 u_3) - \partial_1 p + \Delta u_1 \\
& \partial_t u_2 = -R \partial_1 (u_2 u_1) - R \partial_2 (u_2^2) - R \partial_2 (u_2 u_3) - \partial_2 p + \Delta u_2 \\
& \partial_t u_3 = -R \partial_1 (u_3 u_1) - R \partial_2 (u_3 u_2) - R \partial_2 (u_3^2) - \partial_2 p + \Delta u_3 \\
& \partial_t \rho = - \partial_1 u_1 - \partial_2 u_2 - \partial_3 u_3 \\
& p = \rho/\delta
\end{aligned}
\]

The discretization is as follows:
\[
\begin{aligned}
& u_1^n(i,j,k) = u_1^n(i \Delta x_1, j \Delta x_1, k \Delta x_3, n \Delta t) \\
& u_2^n(i,j,k) = u_2(i \Delta x_1, j \Delta x_1, k \Delta x_3, n \Delta t) \\
& u_3^n(i,j,k) = u_3(i \Delta x_1, j \Delta x_1, k \Delta x_3, n \Delta t) \\
& \rho^n(i,j,k) = \rho(i \Delta x_1, j \Delta x_1, k \Delta x_3, n \Delta t) \\
\end{aligned}
\]

\[
\begin{aligned}
\dfrac{u_1^{n+1}(i,j,k) - u_1^{n-1}(i,j,k)}{2 \Delta t} =& - R \dfrac{(u_1^n(i+1,j,k))^2 - (u_1^n(i-1,j,k))^2}{2 \Delta x_1} \\
& - R \dfrac{u_1^n(i,j+1,k)u_2^n(i,j+1,k) - u_1^n(i,j-1,k)u_2^n(i,j-1,k)}{2 \Delta x_2} \\
& - R \dfrac{u_1^n(i,j,k+1)u_3^n(i,j,k+1) - u_1^n(i,j,k-1)u_3^n(i,j,k-1)}{2 \Delta x_3} \\
& - \dfrac{1}{\delta} \dfrac{\rho^{n+1}(i + 1,j,k) - \rho^n(i-1,j,k)}{2 \Delta x_1} \\
& + \dfrac{u_1^n(i+1,j,k) - 2u_1^n(i,j,k) + u_1^n(i-1,j,k)}{\Delta x_1^2} \\
& + \dfrac{u_1^n(i,j+1,k) - 2u_1^n(i,j,k) + u_1^n(i,j-1,k)}{\Delta x_2^2} \\
& + \dfrac{u_1^n(i,j,k+1) - 2u_1^n(i,j,k) + u_1^n(i,j,k-1)}{\Delta x_3^2}
\end{aligned}
\]

\[
\begin{aligned}
\dfrac{u_2^{n+1}(i,j,k) - u_2^{n-1}(i,j,k)}{2 \Delta t} =& - R \dfrac{u_2^n(i+1,j,k)u_1^n(i+1,j,k) - u_2^n(i-1,j,k)u_1^n(i-1,j,k)}{2 \Delta x_1} \\
& - R \dfrac{(u_2^n(i,j+1,k))^2 - (u_2^n(i,j-1,k))^2}{2 \Delta x_2} \\
& - R \dfrac{u_2^n(i,j,k+1)u_3^n(i,j,k+1) - u_2^n(i,j,k-1)u_3^n(i,j,k-1)}{2 \Delta x_3} \\
& - \dfrac{1}{\delta} \dfrac{p^{n+1}(i,j+1,k) - p^n(i,j-1,k)}{2 \Delta x_2} \\
& + \dfrac{u_2^n(i+1,j,k) - 2u_2^n(i,j,k) + u_2^n(i-1,j,k)}{\Delta x_1^2} \\
& + \dfrac{u_2^n(i,j+1,k) - 2u_2^n(i,j,k) + u_2^n(i,j-1,k)}{\Delta x_2^2} \\
& + \dfrac{u_2^n(i,j,k+1) - 2u_2^n(i,j,k) + u_2^n(i,j,k-1)}{\Delta x_3^2}.
\end{aligned}
\]

\[
\begin{aligned}
\dfrac{u_3^{n+1}(i,j,k) - u_3^{n-1}(i,j,k)}{2 \Delta t} = & - R \dfrac{u_3^n(i+1,j,k)u_1^n(i+1,j,k) - u_3^n(i-1,j,k)u_1^n(i-1,j,k)}{2 \Delta x_1} \\
& -R \dfrac{u_3^n(i,j+1,k)u_2^n(i,j+1,k) - u_3^n(i,j-1,k)u_2^n(i,j-1,k)}{2 \Delta x_2} \\
& - R \dfrac{(u_3^n(i,j,k+1))^2 - (u_3^n(i,j,k-1))^2}{2 \Delta x_3} \\
& - \dfrac{1}{\delta} \dfrac{p^{n+1}(i,j,k+1) - p^n(i,j,k-1)}{2 \Delta x_3} \\
& + \dfrac{u_3^n(i+1,j,k) - 2 u_3^n(i,j,k) + u_3^n(i-1,j,k)}{\Delta x_1^2} \\
& + \dfrac{u_3^n(i,j+1,k) - 2 u_3^n(i,j,k) + u_3^n(i,j-1,k)}{\Delta x_2^2} \\
& + \dfrac{u_3^n(i,j,k+1) - 2 u_3^n(i,j,k) + u_3^n(i,j,k-1)}{\Delta x_3^2}.
\end{aligned}
\]

Therefore
\[
\begin{aligned}
u_1^{n+1}(i,j,k) - u_1^{n-1}(i,j,k) =& - R \dfrac{\Delta t}{\Delta x_1}[(u_1^n(i+1,j,k))^2 - (u_1^n(i-1,j,k))^2] \\
& - R \dfrac{\Delta t}{\Delta x_2}[u_1^n(i,j+1,k)u_2^n(i,j+1,k) - u_1^n(i,j-1,k)u_2^n(i,j-1,k)] \\
& - R \dfrac{\Delta t}{\Delta x_3}[u_1^n(i,j,k+1)u_3^n(i,j,k+1) - u_1^n(i,j,k-1)u_3^n(i,j,k-1)] \\
& - R \dfrac{\Delta t}{\delta \Delta x_1}[\rho^{n+1}(i + 1,j,k) - \rho^n(i-1,j,k)] \\
& + \dfrac{2\Delta t}{\Delta x_1^2}[u_1^n(i+1,j,k) - 2u_1^n(i,j,k) + u_1^n(i-1,j,k)] \\
& + \dfrac{2\Delta t}{\Delta x_1^2} [u_1^n(i,j+1,k) - 2u_1^n(i,j,k) + u_1^n(i,j-1,k)] \\
& + \dfrac{2\Delta t}{\Delta x_1^2}[u_1^n(i,j,k+1) - 2u_1^n(i,j,k) + u_1^n(i,j,k-1)]
\end{aligned}
\]


\[
\begin{aligned}
u_2^{n+1}(i,j,k) - u_2^{n-1}(i,j,k) =& - R \dfrac{\Delta t}{\Delta x_1} [u_2^n(i+1,j,k)u_1^n(i+1,j,k) - u_2^n(i-1,j,k)u_1^n(i-1,j,k)] \\
& - R \dfrac{\Delta t}{\Delta x_2} [(u_2^n(i,j+1,k))^2 - (u_2^n(i,j-1,k))^2] \\
& - R \dfrac{\Delta t}{\Delta x_3} [u_2^n(i,j,k+1)u_3^n(i,j,k+1) - u_2^n(i,j,k-1)u_3^n(i,j,k-1)] \\
& - \dfrac{\Delta t}{\delta \Delta x_2} [p^{n+1}(i,j+1,k) - p^n(i,j-1,k)] \\
& + \dfrac{2\Delta t}{\Delta x_1^2}[u_2^n(i+1,j,k) - 2u_2^n(i,j,k) + u_2^n(i-1,j,k)] \\
& + \dfrac{2\Delta t}{\Delta x_2^2}[u_2^n(i,j+1,k) - 2u_2^n(i,j,k) + u_2^n(i,j-1,k)] \\
& + \dfrac{2\Delta t}{\Delta x_3^2}[u_2^n(i,j,k+1) - 2u_2^n(i,j,k) + u_2^n(i,j,k-1)]
\end{aligned}
\]

\[
\begin{aligned}
u_3^{n+1}(i,j,k) - u_3^{n-1}(i,j,k) = & - R \dfrac{\Delta t}{\Delta x_1} [u_3^n(i+1,j,k)u_1^n(i+1,j,k) - u_3^n(i-1,j,k)u_1^n(i-1,j,k)] \\
& -R \dfrac{\Delta t}{\Delta x_2}[u_3^n(i,j+1,k)u_2^n(i,j+1,k) - u_3^n(i,j-1,k)u_2^n(i,j-1,k)] \\
& - R \dfrac{\Delta t}{\Delta x_3}[(u_3^n(i,j,k+1))^2 - (u_3^n(i,j,k-1))^2] \\
& - \dfrac{\Delta t}{\delta \Delta x_3}[p^{n+1}(i,j,k+1) - p^n(i,j,k-1)] \\
& + \dfrac{2\Delta t}{\Delta x_1^2} [u_3^n(i+1,j,k) - 2 u_3^n(i,j,k) + u_3^n(i-1,j,k)] \\
& + \dfrac{2\Delta t}{\Delta x_2^2} [u_3^n(i,j+1,k) - 2 u_3^n(i,j,k) + u_3^n(i,j-1,k)] \\
& + \dfrac{2\Delta t}{\Delta x_3^2} [u_3^n(i,j,k+1) - 2 u_3^n(i,j,k) + u_3^n(i,j,k-1)]
\end{aligned}
\]


\[
\begin{aligned}
\rho^{n+1}(i,j,k) - \rho^n(i,j,k) =& - \dfrac{\Delta t}{2 \Delta x_1} [u_1^n(i+1,j,k) - u_1^n(i-1,j,k)] \\
& - \dfrac{\Delta t}{2 \Delta x_2} [u_2^n(i+1,j,k) - u_2^n(i-1,j,k)] \\
& - \dfrac{\Delta t}{2 \Delta x_3} [u_3^n(i+1,j,k) - u_3^n(i-1,j,k)]
\end{aligned}
\]


\[
\begin{aligned}
& u_{ji}^{n+1}\left(\dfrac{\Delta x \Delta y}{\Delta t} + \dfrac{2}{Re}\dfrac{\Delta y}{\Delta x} + \dfrac{2}{Re}\dfrac{\Delta x}{\Delta y}\right) - \dfrac{\Delta x \Delta y}{\Delta t}u_{jk}^n + \Delta y \left(p_{j+1,k}^{n+1} - p_{j,k}^{n+1}\right) - \dfrac{1}{Re} \left(u_{j+1,k}^{n+1} + u_{j-1,k}^{n+1} + u_{j,k+1}^{n+1} + u_{j,k-1}^{n+1}\right) \\
& + \Delta y\left[\dfrac{1}{4}(u_{j+1,k}^{n} + u_{j,k}^{n})(u_{j+1,k}^{n+1} + u_{j,k}^{n+1}) - \dfrac{1}{4}(u_{j-1,k}^{n} + u_{j,k}^{n})(u_{j-1,k}^{n+1} + u_{j,k}^{n+1})\right] + \Delta x \left[\dfrac{1}{4}(v_{j,k}^{n} + v_{j+1,k}^{n})(u_{j,k}^{n+1} + u_{j,k+1}^{n+1}) - \dfrac{1}{4}(v_{j,k-1}^{n} + v_{j+1,k-1}^{n})(u_{j,k-1}^{n+1} + u_{j,k}^{n+1})\right] = 0
\end{aligned}
\]

*****************************************

\[
\begin{aligned}
& \left[\dfrac{\Delta x \Delta y}{\Delta t} + \dfrac{2}{Re}\left(\dfrac{\Delta x}{\Delta y} + \dfrac{\Delta y}{\Delta x}\right) + \dfrac{\Delta y}{4}\left(u_{j+1,k}^n - u_{j-1,k}^n\right) + \dfrac{\Delta x}{4}\left(v_{j,k}^n + v_{j+1,k}^n - v_{j,k-1}^n - v_{j+1,k-1}^n\right)\right] u_{ji}^{n+1} \\
& + \left[-\dfrac{1}{Re} + \dfrac{\Delta y}{4}\left(u_{j+1,k}^n + u_{j,k}^n\right)\right] u_{j+1,k}^{n+1} + \left[-\dfrac{1}{Re} - \dfrac{\Delta y}{4}\left(u_{j-1,k}^n + u_{j,k}^n\right)\right] u_{j-1,k}^{n+1} + \left[-\dfrac{1}{Re} + \dfrac{\Delta x}{4}\left(v_{j,k}^n + v_{j+1,k}^n\right)\right] u_{j,k+1}^{n+1} + \left[-\dfrac{1}{Re} - \dfrac{\Delta x}{4}\left(v_{j,k-1}^n + v_{j+1,k-1}^n\right)\right] u_{j,k-1}^{n+1} \\
& - \dfrac{\Delta x \Delta y}{\Delta t} u_{jk}^n + \Delta y \left(p_{j+1,k}^{n+1} - p_{j,k}^{n+1}\right) = 0
\end{aligned}
\]

AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

\[
\begin{aligned}
& \left[\dfrac{\Delta x \Delta y \Delta z}{\Delta t} + \dfrac{2}{Re} \left(\dfrac{\Delta x}{\Delta y} + \dfrac{\Delta y}{\Delta x} + \dfrac{\Delta z}{\Delta x} + \dfrac{\Delta z}{\Delta y} + \dfrac{\Delta x}{\Delta z} + \dfrac{\Delta y}{\Delta z}\right) + \dfrac{\Delta y}{4} \left(u_{j+1,k,l}^n - u_{j-1,k,l}^n\right) \right. \\
& \quad \left. + \dfrac{\Delta x}{4} \left(v_{j,k+1,l}^n - v_{j,k-1,l}^n\right) + \dfrac{\Delta z}{4} \left(w_{j,k,l+1}^n - w_{j,k,l-1}^n\right)\right] u_{j,k,l}^{n+1} \\
& + \left[-\dfrac{1}{Re} + \dfrac{\Delta y}{4} \left(u_{j+1,k,l}^n + u_{j,k,l}^n\right)\right] u_{j+1,k,l}^{n+1} + \left[-\dfrac{1}{Re} - \dfrac{\Delta y}{4} \left(u_{j-1,k,l}^n + u_{j,k,l}^n\right)\right] u_{j-1,k,l}^{n+1} \\
& + \left[-\dfrac{1}{Re} + \dfrac{\Delta x}{4} \left(v_{j,k+1,l}^n + v_{j,k,l}^n\right)\right] u_{j,k+1,l}^{n+1} + \left[-\dfrac{1}{Re} - \dfrac{\Delta x}{4} \left(v_{j,k-1,l}^n + v_{j,k,l}^n\right)\right] u_{j,k-1,l}^{n+1} \\
& + \left[-\dfrac{1}{Re} + \dfrac{\Delta z}{4} \left(w_{j,k,l+1}^n + w_{j,k,l}^n\right)\right] u_{j,k,l+1}^{n+1} + \left[-\dfrac{1}{Re} - \dfrac{\Delta z}{4} \left(w_{j,k,l-1}^n + w_{j,k,l}^n\right)\right] u_{j,k,l-1}^{n+1} \\
& - \dfrac{\Delta x \Delta y \Delta z}{\Delta t} u_{j,k,l}^n + \Delta y \left(p_{j+1,k,l}^{n+1} - p_{j,k,l}^{n+1}\right) = 0.
\end{aligned}
\]
