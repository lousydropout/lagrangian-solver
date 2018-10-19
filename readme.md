# A Lagrangian solver

- *keywords*: Lagrangian, automatic differentiation, backpropagation, ODEs, PDEs
- *total lines of code*:
  - *Ada (2012)*: 7,933
  - *C*: 63

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

Beyond the two, I've also had to write a basic library for dealing with *sparse matrices*. The majority of what I had to write is the basic addition, subtraction, and multiplication of sparse matrices and vectors as well as a wrapper to a sparse linear solver. The sparse linear solver I made use of is called *CXsparse*, which is part of the [SuiteSparse library](http://faculty.cse.tamu.edu/davis/suitesparse.html).
