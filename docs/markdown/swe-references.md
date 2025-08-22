# Bibliography

[◀ Back to Table of Contents](README.md)

## Conclusions

Here lies the shallow-water equations numerical model as described in Kantha and Clayson (2000) with the same numerical scheme. It currently only has the Dirichelet boundary conditions and the gravity wave explicit radiation scheme added of a null-gradient or Flather for the normal velocity. This means that, in the former case, all surface waves bounce back at the boundary and and give rise to a cascade of multiple linear superpositions leading to a path of unavoidable numerical instability. In the latter case, the solution radiates any level perturbation (gravitic waves) at the open boundaries (Orlanski, 1976; Shchepetkin and McWilliams, 2003). 

It is important to note that the geometric considerations of the gaussian bump elevation at the instant of release were crucial in order to estimate matching predictions of energy partitioning and production of TKE. In particular, it was found that the relative production rate of TKE varies with viscosity and $\sigma$ alone, and is independent of the Froude number associated with the gaussian gravity wave.

Further work involves implementing relaxing conditions at the boundaries, variable coriolis force, cyclic boundary conditions. Ultimately, using a sponge layer near the boundaries is considered (Martinsen and Engedahl, 1987; Shchepetkin and McWilliams, 2003; Pietrzak et al., 2002), as well as developing the recent works of Blayo and Debreu (2005) with incoming characteristics and of Lavelle and Thacker (2008) with the PML conditions.

## References

- Arakawa, A. (1966). Computational design for long-term numerical integration of the equations of fluid motion: Two-dimensional incompressible flow. Part I. Journal of Computational Physics, 1(1), 119-143.
- Asselin, R. (1972). Frequency filter for time integrations. Monthly Weather Review, 100(6), 487-490.
- Beckmann, A., & Haidvogel, D. B. (1993). Numerical simulation of flow around a tall isolated seamount. Part I: Problem formulation and model accuracy. Journal of Physical Oceanography, 23(8), 1736-1753.
- Blayo, E., & Debreu, L. (2005). Revisiting open boundary conditions from the point of view of characteristic variables. Ocean modelling, 9(3), 231-252.
- Burchard, H. (2002). Applied turbulence modelling in marine waters. Springer Science & Business Media.
- Courant, R., Friedrichs, K., & Lewy, H. (1959). On the partial difference equations of mathematical physics. IBM Journal of Research and Development, 3(3), 215-234.
- Flather, R. A. (1976). A tidal model of the northwest European continental shelf. Mémoires de la Société Royale des Sciences de Liège, 10(6), 141-164.
- Gill, A. E. (1982). Atmosphere-ocean dynamics (Vol. 30). Academic press.
- Herzfeld, M. (2008). The role of numerical implementation on open boundary behaviour in limited area ocean models. Ocean Modelling, 27(1-2), 18-32.
- Isern-Fontanet, J., Font, J., García-Ladona, E., Emelianov, M., Millot, C., & Taupier-Letage, I. (2004). Spatial structure of anticyclonic eddies in the Algerian basin (Mediterranean Sea) analyzed using the Okubo–Weiss parameter. Deep Sea Research Part II: Topical Studies in Oceanography, 51(25-26), 3009-3028.
- Kantha, L. H., & Clayson, C. A. (2000). Numerical models of oceans and oceanic processes (Vol. 66). Academic press.
- Kundu, P. K., & Cohen, I. M. (2002). Fluid mechanics. Elsevier Academic Press, New York.
- Lavelle, J. W., & Thacker, W. C. (2008). A pretty good sponge: Dealing with open boundaries in limited-area ocean models. Ocean Modelling, 20(3), 270-292.
- Leitao, P. (2003). Modelo de circulação na região do cabo de Sines (in Portuguese). PhD thesis, Technical University of Lisbon.
- Marsaleix, P., Auclair, F., & Estournel, C. (2009). Low-order pressure gradient schemes in sigma coordinate models: The seamount test revisited. Ocean Modelling, 30(2-3), 169-177.
- Martinsen, E. A., & Engedahl, H. (1987). Implementation and testing of a lateral boundary scheme as an open boundary condition in a barotropic ocean model. Coastal engineering, 11(5-6), 603-627.
- Orlanski, I. (1976). A simple boundary condition for unbounded hyperbolic flows. Journal of computational physics, 21(3), 251-269.
- Pedlosky, J. (1987). Geophysical fluid dynamics. Springer, New York.
- Pietrzak, J., Jakobson, J. B., Burchard, H., Vested, H. J., & Petersen, O. (2002). A three-dimensional hydrostatic model for coastal and ocean modelling using a generalised topography following co-ordinate system. Ocean Modelling, 4(2), 173-205.
- Shchepetkin, A. F., & McWilliams, J. C. (2003). A method for computing horizontal pressure-gradient force in an oceanic model with a nonaligned vertical coordinate. Journal of Geophysical Research: Oceans, 108(C3).
- Weiss, J. (1981). The dynamics of enstrophy transfer in two-dimensional hydrodynamics. Technical Report LJI-TN-121, La Jolla Inst., La Jolla, California.
