setwd ("1"); p <- read.csv ("poincare.csv"); # ws = 4.6
## plot (s_dot ~ s, subset (p, pt > 0), cex = 0.1, xlim = c (-1, 1),
##       ylim = c (-4.3, -3), main = "alpha = 10, E = -6.5, (t0, s0) = (pi, -pi)");
plot (s_dot ~ s, subset (p, pt > 0), cex = 0.1, xlim = c (-2, 2),
      ylim = c (-5, 5), main = "alpha = 10, E = -6.5, (t0, s0) = (pi, -pi)");
points (-s_dot ~ s, subset (p, pt < 0), cex = 0.1);
      
setwd ("../2"); p <- read.csv ("poincare.csv"); # ws = 4
points (s_dot ~ s, subset (p, pt > 0), cex = 0.1)
points (-s_dot ~ s, subset (p, pt < 0), cex = 0.1);

setwd ("../3"); p <- read.csv ("poincare.csv"); # ws = 3.5
points (s_dot ~ s, subset (p, pt > 0), cex = 0.1)
points (-s_dot ~ s, subset (p, pt < 0), cex = 0.1);

setwd ("../4"); p <- read.csv ("poincare.csv"); # ws = 3
points (s_dot ~ s, subset (p, pt > 0), cex = 0.1)
points (-s_dot ~ s, subset (p, pt < 0), cex = 0.1)

setwd ("../5"); p <- read.csv ("poincare.csv"); # ws = 2
points (s_dot ~ s, subset (p, pt > 0), cex = 0.1)
points (-s_dot ~ s, subset (p, pt < 0), cex = 0.1);

setwd ("../6"); p <- read.csv ("poincare.csv"); # ws = 1
points (s_dot ~ s, subset (p, pt > 0), cex = 0.1)
points (-s_dot ~ s, subset (p, pt < 0), cex = 0.1);

setwd ("../7"); p <- read.csv ("poincare.csv"); # ws = 2.5
points (s_dot ~ s, subset (p, pt > 0), cex = 0.1)
points (-s_dot ~ s, subset (p, pt < 0), cex = 0.1);

setwd ("../8"); p <- read.csv ("poincare.csv"); # ws = 0
points (s_dot ~ s, subset (p, pt > 0), cex = 0.1)
points (-s_dot ~ s, subset (p, pt < 0), cex = 0.1);

setwd ("../9"); p <- read.csv ("poincare.csv"); # 2.8
points (s_dot ~ s, subset (p, pt > 0), cex = 0.1)
points (-s_dot ~ s, subset (p, pt < 0), cex = 0.1)

setwd ("../10"); p <- read.csv ("poincare.csv"); # -0.5
points (s_dot ~ s, subset (p, pt > 0), cex = 0.1)
points (-s_dot ~ s, subset (p, pt < 0), cex = 0.1);

## setwd ("../11"); p <- read.csv ("poincare.csv"); # 2.9
## points (s_dot ~ s, subset (p, pt > 0), cex = 0.1, col = "red")

setwd ("../12"); p <- read.csv ("poincare.csv"); # -2
points (s_dot ~ s, subset (p, pt > 0), cex = 0.1)
points (-s_dot ~ s, subset (p, pt < 0), cex = 0.1);

setwd ("../13"); p <- read.csv ("poincare.csv"); # 3.1
points (s_dot ~ s, subset (p, pt > 0), cex = 0.1)
points (-s_dot ~ s, subset (p, pt < 0), cex = 0.1)
setwd ("..")
