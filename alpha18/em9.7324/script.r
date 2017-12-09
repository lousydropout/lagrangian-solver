setwd ("1"); p <- read.csv ("poincare.csv");
plot (s_dot ~ s, subset (p, pt > 0), cex = 0.2,
      xlim = c (-0.09, 0.09), ylim = c (-0.4, 0.4),
      main = "alpha = 18, E = -9.7324");
setwd ("../2"); p <- read.csv ("poincare.csv");
points (s_dot ~ s, subset (p, pt > 0), cex = 0.2)
setwd ("../2"); p <- read.csv ("poincare.csv");
points (s_dot ~ s, subset (p, pt > 0), cex = 0.2)
setwd ("../3"); p <- read.csv ("poincare.csv");
points (s_dot ~ s, subset (p, pt > 0), cex = 0.2)
setwd ("../4"); p <- read.csv ("poincare.csv");
points (s_dot ~ s, subset (p, pt > 0), cex = 0.2)
setwd ("../5"); p <- read.csv ("poincare.csv");
points (s_dot ~ s, subset (p, pt > 0), cex = 0.2)
## setwd ("../6"); p <- read.csv ("poincare-1.csv");
## points (s_dot ~ s, subset (p, pt > 0), cex = 0.2)
setwd ("..")
