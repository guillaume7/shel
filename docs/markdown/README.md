# SWAM, a Shallow WAter numerical Model: Scientific Documentation

**Author: Guillaume Riflet**

## Table of Contents

### [Part 1: Model Description](swe-part1.md)
- [Abstract](swe-part1.md#abstract)
- [Keywords](swe-part1.md#keywords)
- [Introduction](swe-part1.md#introduction)
- [Mathematical Model of the Shallow-Water Equations](swe-part1.md#mathematical-model-of-the-shallow-water-equations)
  - [The Mathematical Model](swe-part1.md#the-mathematical-model)
  - [The Mesh](swe-part1.md#the-mesh)
  - [Boundary Conditions](swe-part1.md#boundary-conditions)
    - [Null-flux](swe-part1.md#null-flux)
    - [No-slip](swe-part1.md#no-slip)
    - [Radiative Boundary Conditions](swe-part1.md#radiative-boundary-conditions)
  - [The Numerical Scheme](swe-part1.md#the-numerical-scheme)

### [Part 2: Model Validation](swe-part2.md)
- [Validation](swe-part2.md#validation)
  - [Gaussian Bell-Shaped Geometry](swe-part2.md#gaussian-bell-shaped-geometry)
  - [Energy](swe-part2.md#energy)
  - [Geometric and Similitude Considerations](swe-part2.md#geometric-and-similitude-considerations)
  - [Basic Results](swe-part2.md#basic-results)

### [Part 3: Advanced Analysis and Conclusions](swe-part3.md)
- [Energy Decay Study](swe-part3.md#energy-decay-study)
- [Radiation Boundary Condition](swe-part3.md#radiation-boundary-condition)
- [Geostrophic Equilibrium](swe-part3.md#geostrophic-equilibrium)
- [Applying the Okubo-Weiss Scalar to Assess the Open-Boundary Condition](swe-part3.md#applying-the-okubo-weiss-scalar-to-assess-the-open-boundary-condition)
- [Conclusions](swe-part3.md#conclusions)
- [References](swe-part3.md#references)

## Document Overview

This documentation is split into three parts for better readability:

1. **Part 1** covers the mathematical and numerical fundamentals of the model, including the equations, mesh design, and boundary conditions.
2. **Part 2** focuses on model validation using a Gaussian bump test case, analyzing conservation properties and energy behavior.
3. **Part 3** explores advanced topics like energy decay, radiation boundary conditions, geostrophic equilibrium, and the application of the Okubo-Weiss scalar.

All equations are rendered using standard Markdown equation syntax compatible with both GitHub and VS Code viewers.
