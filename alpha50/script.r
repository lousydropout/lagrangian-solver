setwd ("1"); p <- read.csv ("poincare.csv");
plot (s_dot ~ s, subset (p, pt > 0), cex = 0.2, xlim = c (-1, 1),
      ylim = c (-5, 5), main = "alpha = 50, E = -28");
## points (s_dot ~ s, t_dot < 0), cex = 0.2, col = "red");
setwd ("../2"); p <- read.csv ("poincare.csv");
points (s_dot ~ s, subset (p, pt > 0), cex = 0.2)
## points (s_dot ~ s, subset (p, t_dot < 0), cex = 0.2, col = "red");
setwd ("../2"); p <- read.csv ("poincare.csv");
points (s_dot ~ s, subset (p, pt > 0), cex = 0.2)
## points (s_dot ~ s, subset (p, t_dot < 0), cex = 0.2, col = "red");
setwd ("../3"); p <- read.csv ("poincare.csv");
points (s_dot ~ s, subset (p, pt > 0), cex = 0.2)
## points (s_dot ~ s, subset (p, t_dot < 0), cex = 0.2, col = "red");
setwd ("../4"); p <- read.csv ("poincare.csv");
points (s_dot ~ s, subset (p, pt > 0), cex = 0.2)
## points (s_dot ~ s, subset (p, t_dot < 0), cex = 0.2, col = "red");
setwd ("../5"); p <- read.csv ("poincare.csv");
points (s_dot ~ s, subset (p, pt > 0), cex = 0.2)
## points (s_dot ~ s, subset (p, t_dot < 0), cex = 0.2, col = "red");
setwd ("../6"); p <- read.csv ("poincare.csv");
points (s_dot ~ s, subset (p, pt > 0), cex = 0.2)
## points (s_dot ~ s, subset (p, t_dot < 0), cex = 0.2, col = "red");
setwd ("../7"); p <- read.csv ("poincare.csv");
points (s_dot ~ s, subset (p, pt > 0), cex = 0.2)
## points (s_dot ~ s, subset (p, t_dot < 0), cex = 0.2, col = "red");
setwd ("../8"); p <- read.csv ("poincare.csv");
points (s_dot ~ s, subset (p, pt > 0), cex = 0.2)
## points (s_dot ~ s, subset (p, t_dot < 0), cex = 0.2, col = "red");
setwd ("../9"); p <- read.csv ("poincare.csv");
points (s_dot ~ s, subset (p, pt > 0), cex = 0.2)
## points (s_dot ~ s, subset (p, t_dot < 0), cex = 0.2, col = "red");
setwd ("..")
