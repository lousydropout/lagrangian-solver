setwd ("1"); p <- read.csv ("poincare.csv");
plot (s_dot ~ s, subset (p, pt > 0), cex = 0.2,
      xlim = c (-0.1, 0.1), ylim = c (-0.5, 0.5), main = "alpha = 18, E = -9.74");
setwd ("../2"); p <- read.csv ("poincare-1.csv");
points (s_dot ~ s, subset (p, pt > 0), cex = 0.2)
setwd ("../2"); p <- read.csv ("poincare.csv");
points (s_dot ~ s, subset (p, pt > 0), cex = 0.2)
setwd ("../3"); p <- read.csv ("poincare.csv");
points (s_dot ~ s, subset (p, pt > 0), cex = 0.2)
setwd ("../4"); p <- read.csv ("poincare.csv");
points (s_dot ~ s, subset (p, pt > 0), cex = 0.2)
setwd ("../5"); p <- read.csv ("poincare.csv");
points (s_dot ~ s, subset (p, pt > 0), cex = 0.2)
setwd ("../6"); p <- read.csv ("poincare.csv");
points (s_dot ~ s, subset (p, pt > 0), cex = 0.2)
setwd ("../7"); p <- read.csv ("poincare.csv");
points (s_dot ~ s, subset (p, pt > 0), cex = 0.2)
setwd ("../8"); p <- read.csv ("poincare.csv");
points (s_dot ~ s, subset (p, pt > 0), cex = 0.2)
setwd ("..")
