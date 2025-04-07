import random
import numpy as np
import matplotlib.pyplot as plt
from enum import Enum

# define the DIRECTION enum to simulate random walk directions
class DIRECTION(Enum):
    LEFT = 0
    RIGHT = 1
    UP = 2
    DOWN = 3
    DIRECTION_COUNT = 4

# generate a random walk direction using a uniform distribution
def generateRandomWalkDirection():
    return DIRECTION(random.randint(0, DIRECTION.DIRECTION_COUNT.value - 1))

# function phi defines the boundary conditions for the solution
def phi(x, y, L):
    if x == 0.0:     # boundary on the left edge
        return 1.0 - y     # boundary on the left edge
    elif y == 0.0:
        return 1.0 - x    # boundary on the bottom edge
    elif x == L or y == L:
        return 0.0    # boundary on the right and top edges
    else:
        return -1.0    # EXCEPTION CASE, JUST TO COMPLETE THE CODE AND AVOID WARNING SIGNS

# Check if a point (i, j) is on the boundary of the grid
def isBoundaryPoint(i, j, N):
    if i == 0 or i == N or j == 0 or j == N:
        return True
    else:
        return False

# initialize the boundary values of U using the phi function
def initializeBoundaryNode(U, N, L, dx, dy):
    for j in range(0, N + 1):
        U[0][j] = phi(0 * dx, j * dy, L)    # bottom boundary
        U[N][j] = phi(N * dx, j * dy, L)    # top boundary
        U[j][0] = phi(j * dx, 0 * dy, L)    # left boundary
        U[j][N] = phi(j * dx, N * dy, L)    # right boundary

# count the number of hits at each boundary point during random walks
# @return NOTHING
def countBoundaryPointHit(i, j, C, N):
    if i == 0 and 1 <= j <= N: # bottom boundary
       C[0][j - 1] += 1
    elif i == N and 1 <= j <= N: # top boundary
       C[1][j - 1] += 1
    elif j == 0 and 1 <= i <= N: # left boundary
       C[2][i - 1] += 1
    elif j == N and 1 <= i <= N: # right boundary
       C[3][i - 1] += 1

# solve the Laplace equation using the Monte Carlo random walk method
# @input N, CUSTOMIZABLE, number of rows and columns for the mesh
# @input K, CUSTOMIZABLE, controlled number of random walks
# @return the matrix U[i, j]
def solveLaplaceEquationMarkov(N, K, L):
    U = [[0.0 for _ in range(N + 1)] for _ in range(N + 1)] # initialize the solution matrix U(i, j)
    initializeBoundaryNode(U, N, L, L/float(N), L/float(N))    # initialize the boundary nodes
    C = [[0 for _ in range(N - 1)] for _ in range(4)] # initialize the counters for boundary point hits
    for i in range(1, N):    # loop through the grid points, not including boundary points
        for j in range(1, N):
            # reset the counter for the current cell
            C = [[0 for _ in range(N - 1)] for _ in range(4)]
            for _ in range(K):    # perform K random walks
                currentI, currentJ = i, j
                    # simulate the random walk until hitting a boundary
                while not isBoundaryPoint(currentI, currentJ, N):
                    nextMove = generateRandomWalkDirection()
                    if nextMove == DIRECTION.LEFT and currentI > 0:
                        currentI -= 1
                    elif nextMove == DIRECTION.RIGHT and currentI < N:
                        currentI += 1
                    elif nextMove == DIRECTION.UP and currentJ > 0:
                        currentJ -= 1
                    elif nextMove == DIRECTION.DOWN and currentJ < N:
                        currentJ += 1
                countBoundaryPointHit(currentI, currentJ, C, N)
            # WEIGHTED-MEAN of U[i][j] based on contributions from the boundary points
            for l in range(N - 1):
                U[i][j] += C[0][l] * U[0][l + 1]    # bottom boundary
                U[i][j] += C[1][l] * U[N][l + 1]    # top boundary
                U[i][j] += C[2][l] * U[l + 1][0]    # left boundary
                U[i][j] += C[3][l] * U[l + 1][N]    # right boundary
            U[i][j] /= float(K)    # calculate the probability p[i,j](B)
    # generate the x and y values for the graphic grid
    x_vals = np.linspace(0, l, N + 1)
    y_vals = np.linspace(0, l, N + 1)

    # create a meshgrid for plotting
    X, Y = np.meshgrid(x_vals, y_vals)

    # convert U into a NumPy array for contour plotting
    solution = np.array(U)

    # plot the solution using a filled contour plot
    plt.figure(figsize=(8, 6))
    plt.contourf(X, Y, solution, levels=50, cmap="coolwarm")    # use coolwarm colormap
    plt.colorbar(label="Value of the function")    # add a colorbar for reference
    plt.title("Approximated laplace Equation Solution via Monte Carlo")    # Plot title
    plt.xlabel("x")    # label for x-axis
    plt.ylabel("y")    # label for y-axis
    plt.show()    # display the plot

# solve the laplace equation
solveLaplaceEquationMarkov(10, 1000, 1.0)