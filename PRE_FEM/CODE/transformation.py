import sympy as sp

# Define variables and functions
x, y, delta_x1, delta_x2, delta_y1 = sp.symbols('x y delta_x1 delta_x2 delta_y1')
k1, beta1, beta = sp.symbols('k1 beta1 beta')
u = sp.Function('u')

# Original equation
original_eq = (
    2 * (u(x + delta_x2, y) - u(x, y)) / (k1 * (1 + k1) * beta1**2) +
    2 * (u(x - delta_x1, y) - u(x, y)) / ((1 + k1) * beta1**2) -
    2 * beta * (u(x, y + delta_y1) - u(x, y)) / ((1 + k1) * beta1) -
    2 * beta * (
        u(x + delta_x2, y + delta_y1) - u(x, y + delta_y1) -
        u(x + delta_x2, y) + u(x, y)
    ) / ((1 + k1) * beta1) +
    ((k1 / (1 + k1)) * beta**2 + 1) *
    (u(x, y + delta_y1) + u(x, y - delta_y1) - 2 * u(x, y))
)

# Transformed equation
transformed_eq = (
    -2 * (1 + k1 + k1 * beta1**2 * (1 + k1 + k1 * beta**2)) * u(x, y) +
    2 * (1 + k1 * beta1 * beta) * u(x + delta_x2, y) +
    2 * k1 * u(x - delta_x1, y) +
    k1 * beta1**2 * (1 + k1 + k1 * beta**2) * u(x, y + delta_y1) +
    k1 * beta1**2 * (1 + k1 + k1 * beta**2) * u(x, y - delta_y1) -
    2 * k1 * beta1 * beta * u(x + delta_x2, y + delta_y1)
)

# Simplify the original equation
simplified_original = sp.simplify(original_eq)

# Simplify the transformed equation
simplified_transformed = sp.simplify(transformed_eq)

# Check equality
is_equal = sp.simplify(simplified_original - simplified_transformed) == 0

# Display results
print("Simplified Original Equation:")
print(simplified_original)
print("\nSimplified Transformed Equation:")
print(simplified_transformed)
print("\nAre the equations equivalent?", is_equal)