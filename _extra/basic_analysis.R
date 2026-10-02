# basic_analysis.R
# Generate some random data
set.seed(42)
x <- rnorm(100)
y <- 2 * x + rnorm(100)

# Perform a basic linear regression
model <- lm(y ~ x)

# Print the summary of the model
print(summary(model))
