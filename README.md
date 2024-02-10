# A Specialized Semismooth Newton Method for Kernel-Based Optimal Transport

These codes provide implementations of solvers for solving kernel-based optimal transport problems using a specialized semi-smooth Newton method. 

# About

Kernel-based optimal transport (OT) estimators offer an alternative, functional estimation procedure to address OT problems from samples. Recent works suggest that these estimators are more statistically efficient than
plug-in (linear programming-based) OT estimators when comparing probability measures in high-dimensions (Vacher et al., 2021). Unfortunately, that statistical benefit comes at a very steep computational price: because
their computation relies on the short-step interior-point method (SSIPM), which comes with a large iteration count in practice, these estimators quickly become intractable w.r.t. sample size n. 

To scale these estimators to larger n, we propose a nonsmooth fixed-point model for the kernel-based OT problem, and show that it can be efficiently solved via a specialized semismooth Newton (SSN) method: We show, exploring the problem’s structure, that the per-iteration cost of performing one SSN step can be significantly reduced in practice. We prove that our SSN method achieves a global convergence rate and a local quadratic conver- gence rate under standard regularity conditions. We show substantial speedups over SSIPM on both synthetic and real datasets.

# Codes

The MATLAB Implementations on Synthetic Data are provided.  

# References

T. Lin, M. Cuturi and M. I. Jordan. A Specialized Semismooth Newton Method for Kernel-Based Optimal Transport. AISTATS'2024. 
