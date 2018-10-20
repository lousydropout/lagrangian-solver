# A Lagrangian solver

- *keywords*: Lagrangian, automatic differentiation, backpropagation, ODEs, PDEs
- *total lines of code*:
  - *Ada (2012)*: 7,933
  - *C*: 63

  1. [Introduction](#introduction)
  2. [Techniques used in this project](#techniques-used-in-this-project)
  3. [How to use this library](#how-to-use-this-library)
  4. [Future](#future)

## Introduction

This is the git repository for a new numerical solver I designed and developed that I simply call a *Lagrangian solver*. The method is derived and explained in a lot more detail in my thesis, which can be found in this directory.

A **Lagrangian** can be thought of as something even more *fundamental* to *Netwon's equations of motion*. Whereas Newton's equations of motion is a system of *ND* 2nd-order ordinary differential equations (ODEs), where *N* is the number of particles in your system and *D* is the number of spatial dimensions, a Lagrangian is a *single* expression from which the *ND* equations can be derived.

In short, given an expression for the Lagrangian, we can derive all of the Newton's equations of motion! If you think of science, and physics in particular, as being on a quest for the minimum number of concepts needed to explain everything, you can understand why the Lagrangian is often thought of as being more fundamental the Netwon's equations of motion.

Even though the Lagrangian is simpler, more elegant, and the thing physicists start from (as opposed to Newton's equations of motion), when we have to calculate things, we often have to derive Newton's equations of motion. Not only is this step tedious and error-prone (due to tediousness), it is, in some sense, the antithesis of the progress made in physics! And that is where my *Lagrangian solver* comes in!

By making use of **automatic differentiation** and the **collocation method**, we can avoid doing the aforementioned progress-reversing step ourselves!

## Techniques used in this project

The main conceptual techniques I make use of in this project are:
- the forward mode of **automatic differentiation** and
- the **collocation method**.

The *forward mode of automatic differentiation* is a way of getting the derivatives (gradients and Hessians in our case) of a function without *explicitly* calculating them. It is a variant of automatic differentiation. Another variant is the *reverse mode of automatic differentiation*, which is perhaps better known in *machine learning* as **backpropagation**.

*Collocation method* is a numerical method for solving differential equations. Strangely, while it is *somewhat* popular in solving *partial differential equations (PDEs)*, I have yet to see it be used for *ordinary differential equations (ODEs)*. And so, far as I'm aware, this is the first use of the collocation method for solving ODEs.

Beyond the two, I've also had to write a basic library for dealing with *sparse matrices*. The majority of what I had to write is the basic addition, subtraction, and multiplication of sparse matrices and vectors as well as a wrapper to a sparse linear solver. The sparse linear solver I made use of is called *CXSparse*, which is part of the [SuiteSparse library](http://faculty.cse.tamu.edu/davis/suitesparse.html).


## How to use this library

First, requirements:
- [Ada 2012 compiler](https://www.adacore.com/download)
- [A C compiler](https://gcc.gnu.org/)
- [CXSparse library (part of the SuiteSparse library)](http://faculty.cse.tamu.edu/davis/suitesparse.html)

For now, the three examples available are:
- a simple pendulum problem,
- a chaotic Henon-Heiles problem, and
- a steel balls problem that's hard to explain solely in words (so I won't).

The driver files for these programs are contained in the directory *drivers* (along with many other driver files that were used for testing and exploration purposes), but to run these programs, follow the steps below.
1. Open a terminal on the parent directory (the one containing the file *main.gpr*).
2. Type the command below to compile all 3 examples listed above.
```bash
    gprbuild -P main.gpr
```
3. Then to run the files, use the commands below. Note that while the bin file *henon_heiles* does not take a file input, the other two do.

#### Chaotic Henon-Heiles:
*Note: This program will take a few minutes.*
```bash
    ./henon_heiles
```
This outputs *out.csv*. The following is the *R* code I use to import and plot the data:
```R
  data <- read.csv ("out.csv")
  plot (q2_dot ~ q2, data, pch = 10, cex = 0.2)
```


#### Pendulum:
```bash
    ./pendulum < param.txt
```
This outputs *pendulum.csv* which you can use, for example, R or python to visualize.

#### Steel Balls:
``` bash
    ./steel_balls < param-steel_balls.txt
```
This will output two files: *sb.csv* and *sb.xyz*.

- *sb.xyz* can be read into *Ovito*, a visualization program mainly meant for molecular dynamics. You can load it into *Ovito* via *Load file (Ctrl+I)*. Then, in the lower-right hand corner, inside the *XYZ* tab, click on *File contains time steps* (else it won't know to look for the particles' positions at other timesteps in the same file). Lastly, if the particles look like they're overlapping or not touching each other, confirm that the particle radius is set to *1*. Then click on the *Play* buttom near the bottom to start the animation.
- *sb.csv* contains more detailed information, but is not a file format readable by *Ovito*.


## Future

Many improvements and extension remain to be done and I'm not sure if I'll ever get to them. Nevertheless, let me note just two here (one of each).
1. The automatic differentiation here is implemented by overloading the various operators such as *+*, *-*, *\**, */*, *and* (for tensor multiplication), and *or* (for direct summation). I don't know what I'm doing wrong, but this seems to be creating a lot of overhead, and the program is a lot slower than I expected.
2. Much like in regular molecular dynamics, by implementing a octree algorithm, it is possible to reduce the number of terms in the Lagrangian and get the same speedup as in molecular dynamics. (Also, I'll need to review the formulations of the various isostats, and implement them, before this will act as a proper molecular dynamics program.)
