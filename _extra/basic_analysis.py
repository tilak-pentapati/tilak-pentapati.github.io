import numpy as np

# Generate some random data
np.random.seed(42)
x = np.random.randn(100)
y = 2 * x + np.random.randn(100)

# Perform a basic linear regression
coefficients = np.polyfit(x, y, 1)

print("Linear Regression Results:")
print(f"Slope (Coefficient): {coefficients[0]:.4f}")
print(f"Intercept: {coefficients[1]:.4f}")
