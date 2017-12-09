setwd ("1"); p <- read.csv ("poincare.csv");
plot (s_dot ~ s, p, cex = 0.2, xlim = c (-0.3, 0.3),
      ylim = c (-1.5, 1.5), main = "alpha = 20, E = -11");
points (s_dot ~ s, subset (p, t_dot < 0), cex = 0.2, col = "red");
setwd ("../2"); p <- read.csv ("poincare.csv");
points (s_dot ~ s, p, cex = 0.2)
points (s_dot ~ s, subset (p, t_dot < 0), cex = 0.2, col = "red");
setwd ("../2"); p <- read.csv ("poincare.csv");
points (s_dot ~ s, p, cex = 0.2)
points (s_dot ~ s, subset (p, t_dot < 0), cex = 0.2, col = "red");
setwd ("../3"); p <- read.csv ("poincare.csv");
points (s_dot ~ s, p, cex = 0.2)
points (s_dot ~ s, subset (p, t_dot < 0), cex = 0.2, col = "red");
setwd ("../4"); p <- read.csv ("poincare.csv");
points (s_dot ~ s, p, cex = 0.2)
points (s_dot ~ s, subset (p, t_dot < 0), cex = 0.2, col = "red");
setwd ("../5"); p <- read.csv ("poincare.csv");
points (s_dot ~ s, p, cex = 0.2)
points (s_dot ~ s, subset (p, t_dot < 0), cex = 0.2, col = "red");
setwd ("..")
