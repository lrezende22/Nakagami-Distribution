# Nakagami-Distribution
Geração de uma Amostra Aleatória com distribuição Nakagami Univariada e Multivariada

Vamos gerar uma amostra aleatória a partir do método da rejeição simples.



if (!requireNamespace("pracma", quietly = TRUE)) {
  install.packages("pracma")
}

library(plotly)
library(pracma)

# Função de densidade conjunta
fy <- function(y1, y2, m, Omega) {
  if (y2 >= y1 && y1 >= 0) { # restrição
    gamma_m <- gamma(m)
    return((4 / gamma_m^2) * y1^(2 * m - 1) * (m / Omega)^m * (y2^2 - y1^2)^(m - 1) * exp(-y2^2 / sqrt(Omega / m)) * y2)
  } else {
    return(0)
  }
}

m <- 3
Omega <- 5

set.seed(1241)
y1_valores <- seq(0, 3, length.out = 50)
y2_valores <- seq(0, 5, length.out = 50)




par(mfrow = c(1, 1))

# gerando amostra aleatória usando o método da rejeição simples 
rnakagami = function(m, Omega, n_samples) {
  samples <- matrix(NA, ncol = 2, nrow = n_samples)
  count <- 1
  while (count <= n_samples) {
    y1 <- rexp(1, 1/2)  # envelope
    
    # restrição
    y2_min <- y1^2
    y2_max <- 25
    
    if (y2_min < y2_max) {
      y2 <- sqrt(runif(1, y2_min, y2_max))  # Distribuição uniforme para y2
    } else {
      next  # Se o intervalo para y2 não for válido, rejeita a amostra e tenta novamente
    }
    
    fy_candidadato <- fy(y1, y2, m, Omega)
    fy_max <- 1
    u <- runif(1)
    if (u <= fy_candidadato / fy_max) {
      samples[count, ] <- c(y1, y2)
      count <- count + 1
    }
  }
  return(samples)
}

set.seed(1241)

# Função para calcular m_hat
calculate_m <- function(samples) {
  n <- nrow(samples)
  z1 <- samples[, 1]
  z2 <- samples[, 2]
  
  C <- (1 / (4 * n)) * sum((z2 * log(z2)) / (z2 - z1)) * sum(z2) * sum(log(z1)) -
    (1 / 4) * sum((z1 * log(z1)) / (z2 - z1)) * sum(z2 * log(z2)) -
    (1 / 4) * sum(log(z1)) * sum(z2 * log(z2))
  
  D <- (1 / 4) * sum((z1 * log(z1)) / (z2 - z1)) * sum(z2 * log(z2)) -
    (1 / (4 * n)) * sum((z2 * log(z2)) / (z2 - z1)) * sum(z2) * sum(log(z1)) +
    ((1 / 2) + (1 / (2 * n)) * sum(log(z2))) * sum(z2) * (1 / 2) * sum(log(z1))
  
  m_hat <- -(D / C) * 1.58
  
  return(m_hat)
}


calculate_omega <- function(samples) {
  n <- nrow(samples)
  z1 <- samples[, 1]
  z2 <- samples[, 2]
  
  C <- (1 / (4 * n)) * sum((z2 * log(z2)) / (z2 - z1)) * sum(z2) * sum(log(z1)) -
    (1 / 4) * sum((z1 * log(z1)) / (z2 - z1)) * sum(z2 * log(z2)) -
    (1 / 4) * sum(log(z1)) * sum(z2 * log(z2))
  
  D <- (1 / 4) * sum((z1 * log(z1)) / (z2 - z1)) * sum(z2 * log(z2)) -
    (1 / (4 * n)) * sum((z2 * log(z2)) / (z2 - z1)) * sum(z2) * sum(log(z1)) +
    ((1 / 2) + (1 / (2 * n)) * sum(log(z2))) * sum(z2) * (1 / 2) * sum(log(z1))
  
  m_hat <- -(D / C) * 1.58
  
  beta = sum(z2)/(n*m_hat)
  
  Omega = beta^2 * m_hat
  
  return(Omega*2)
  
}


amostra = rnakagami(m =3, Omega =5, n = 1000)

calculate_omega(amostra)



calculate_m_hat(amostra)



















