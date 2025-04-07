import sympy as sp

# Define variables and functions
x, y, delta_x2, delta_x3, delta_y3 = sp.symbols('x y delta_x2 delta_x3 delta_y3')
k2, beta3, beta = sp.symbols('k2 beta3 beta')
u = sp.Function('u')

# Original equation
original_eq = (
    2 * (u(x + delta_x3, y) - u(x, y)) / (beta3 * (beta3 / k2 + beta3)) +
    2 * (u(x - delta_x2, y) - u(x, y)) / ((beta3 / k2) * (beta3 / k2 + beta3)) -
    2 * beta * (u(x, y + delta_y3) - u(x, y)) / (beta3 / k2 + beta3) -
    2 * beta * (
        u(x - delta_x2, y + delta_y3) - u(x, y + delta_y3) -
        u(x - delta_x2, y) + u(x, y)
    ) / (beta3 / k2 + beta3) +
    ((beta**2 / (1 + k2)) + 1) *
    (u(x, y + delta_y3) + u(x, y - delta_y3) - 2 * u(x, y))
)

# Transformed equation
transformed_eq = (
    -2 * (k2 + k2**2 + beta3**2 * (1 + k2 + beta**2)) * u(x, y) +
    2 * k2 * u(x + delta_x3, y) +
    2 * k2 * (k2 + beta * beta3) * u(x - delta_x2, y) +
    beta3**2 * (1 + k2 + beta**2) * u(x, y + delta_y3) +
    beta3**2 * (1 + k2 + beta**2) * u(x, y - delta_y3)
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