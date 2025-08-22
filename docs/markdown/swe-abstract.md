# SHEL, a SHallow water equations modEL: Scientific Documentation v0.9

**Author: Guillaume Riflet**

[◀ Back to Table of Contents](README.md)

## Abstract

The SHallow water equations modEL (SHEL) is a development environment allowing the rapid prototyping of new finite-difference numerical schemes. It is suited for undergrad and grad-students who wish to learn how finite-difference numerical schemes are actually implemented and/or wish to implement their own. It is also suited to students (or professors) who simply wish to visualize some simple scenarios of shallow-water flows (for their students). 

SHEL is developed in Matlab and equipped with a Matlab GUI for easy loading, running and visualization of case-studies (it exports in png and eps formats and does avi movies). The model comes with a series of pre-configured test-cases. New test-cases can easily be implemented, saved and shared with peers. The program is built so that other developers can replace fairly easily the built-in numerical schemes with new numerical schemes and, eventually, contribute to the available list of numerical schemes for SHEL (finite-difference-based only). 

By default, the model comes with the shallow water equations discretized in an Arakawa C-grid with variable bottom, free-lid, land-mask, and a leapfrog and central differences scheme combined with simple Asselin-Roberts filtering, as presented in Kantha and Clayson (2000) and further elaborated in this document. Dirichelet, Neummann and Sommerfeld type conditions are implemented at the boundaries. Simple tests were performed with a gaussian level elevation where the conservation of volume, momentum, mechanical energy and vorticity were analyzed.

## Keywords

- Shallow-waters equations
- Open boundary condition
- Okubo-Weiss scalar
- Numerical model
