A <- matrix(c(3, -4, 10, 5), 2, 2, byrow=TRUE)
A
t(A)
res<-eigen(A)
values <- res$values
vectors<- res$vectors 

# Compute dot product
dot_product <- round(vectors[1,] %*% vectors[2,],8)

# Print result
print(dot_product)


# Create a symmetric matrix
A <- matrix(c(1, 2, 1, 3), nrow = 2)
A

# Get eigenvectors
eig <- eigen(A)
vectors <- eig$vectors  # columns are eigenvectors

# Plot setup
plot(NA, xlim = c(-2, 2), ylim = c(-2, 2), asp = 1, xlab = "X", ylab = "Y", main = "Orthogonal Eigenvectors")
abline(h = 0, v = 0, col = "gray")

# Draw the eigenvectors as arrows from the origin
arrows(0, 0, vectors[1, 1], vectors[2, 1], col = "blue", lwd = 2)
arrows(0, 0, vectors[1, 2], vectors[2, 2], col = "red", lwd = 2)

# Add legend
legend("topright", legend = c("v1", "v2"), col = c("blue", "red"), lwd = 2)
