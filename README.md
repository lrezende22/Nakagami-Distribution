----

title: Nakagami-Distribution
output: htlm_document

---- 
## Geração de uma Amostra Aleatória com distribuição Nakagami Univariada e Multivariada

#Vamos gerar uma amostra aleatória a partir do método da rejeição simples.

```{r}
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

# Gerando amostra aleatória usando o método da rejeição simples 
rnakagami <- function(m, Omega, n_samples) {
  samples <- matrix(NA, ncol = 2, nrow = n_samples)
  count <- 1
  while (count <= n_samples) {
    y1 <- rexp(1, 1/2)  
    
    # Restrição
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

amostra = rnakagami(m=3, Omega = 4, n = 100)

head(amostra)

```
