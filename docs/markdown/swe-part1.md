# SHEL, a SHallow water equations modEL: Scientific Documentation v0.9

**Author: Guillaume Riflet**

[◀ Back to Table of Contents](README.md)

## Introduction

The SHallow water equations modEL (SHEL) is a development environment for implementing and testing finite-difference numerical schemes for shallow water equations. It provides a platform for educational purposes, allowing students and researchers to understand how numerical discretization affects the simulation of coastal and oceanic processes.

This document presents the mathematical formulation, numerical implementation, and validation of the model. It covers the core equations, grid structure, boundary conditions, and the numerical schemes used to solve the shallow water equations. The model's performance is assessed through various test cases focused on conservation properties and physical consistency.

## Mathematical Model of the Shallow-Water Equations

### The Mathematical Model

**Fig. 1:** ![System depicted by the mathematical model. The reference level is indicated by a dash-dotted line. $\eta$ is the water elevation from the reference level, $d$ is the depth from the reference level and $H$ is the total depth. The forces acting on the system are illustrated by vector arrows. $g$ is the gravitational acceleration, $\tau_w$ is the wind stress, $\tau_b$ is the bottom stress and $\Omega \times \textbf{v}$ is the Coriolis acceleration.](figs/swe-system-illustration-xana.svg)

The shallow waters equations (SWE) describe the 2D barotropic
motion of water masses. The system and its forcings are
illustrated in Fig. 1. The SWE are widely
described throughout the literature; for example, they are given
in Kantha and Clayson (2000) as:

$$
\begin{cases}
 \frac{\partial Hu}{\partial t} + \frac{\partial \left(Huu\right)}{\partial x} + \frac{\partial \left(Huv\right)}{\partial y} - fHv = \nu \left( \frac{\partial}{\partial x}\left( H \frac{\partial u}{\partial x }\right) + \frac{\partial}{\partial y}\left(H \frac{\partial u}{\partial y }\right) \right) \\
 \quad - gH\frac{\partial \eta}{\partial x} + \frac{\tau _u^w}{\rho _0} - \frac{\tau _u^b}{\rho _0} \\
 \frac{\partial Hv}{\partial t} + \frac{\partial \left(Hvu\right)}{\partial x} + \frac{\partial \left(Hvv\right)}{\partial y} + fHu = \nu \left( \frac{\partial}{\partial x}\left( H \frac{\partial v}{\partial x }\right) + \frac{\partial}{\partial y}\left(H \frac{\partial v}{\partial y }\right) \right) \\
 \quad - gH\frac{\partial \eta}{\partial y} + \frac{\tau _v^w}{\rho _0} - \frac{\tau _v^b}{\rho _0} \\
 \frac{\partial \eta}{\partial t} + \frac{\partial Hu}{\partial x} + \frac{\partial Hv}{\partial y} = 0
\end{cases}
$$

where $H$ is the depth from the surface level to the bottom, $u$
and $v$ are the zonal and meridional components of velocity, $x$,
$y$ and $z$ are the zonal, meridional and depth positions
respectively. $f=1.01\times 10^{-4}$ rad/s is the Coriolis
frequency at 42° of latitude, $\nu$ is the horizontal turbulent
viscosity, $g=9.8$ m$^2$/s is the gravity acceleration,
$\rho_0=1033$ kg/m$^3$ is the water mean density and $\eta$ is the
water level relative to rest. $\tau _u^b$ is the bottom stress
zonal component, $\tau _u^w$ is the wind stress zonal component.
The bottom stress (Pietrzak et al., 2002) is given by:

$$
\tau _u^b = \rho _0 C_D u_b \sqrt {u_b^2 + v_b^2 }
$$

where $C_D$ is the bottom drag coefficient and $u_b$ and $v_b$ are
the zonal and meridional velocity bottom velocity components. The
bottom drag coefficient (Leitao, 2003) is given by:

$$
C_D = \left(k/\ln \left( \frac{z_D + z_0}{z_0} \right)\right)^2
$$

where $z_D$ is the distance to the bottom, $z_0 =0.002$ m is a
typical roughness length (Leitao, 2003) and the Von
Karman constant (Leitao, 2003) is set to:

$$
k = 0.4.
$$

The wind stress (Pietrzak et al., 2002) is given by:

$$
\tau _u^w = \rho _a C_a u_{10} \sqrt {u_{10}^2 + v_{10}^2 }
$$

where $\rho _a = 1.25$ kg/m³ is the air density, $C_a$ is an air
drag coefficient whose values can be found in Pietrzak et al. (2002) and $u_{10}$ and $v_{10}$ is the air speed at
10m height.

### The Mesh

**Fig. 2:** ![Arakawa C staggered grid patterns. From left to right: the T-cell, where $\eta$ and $H$ are evaluated at the centres, and $u$ and $v$ are evaluated at the eastern, western faces and southern, northern faces respectively. The U-cell where $u$ is evaluated at the centre, $\eta$ and $H$ are evaluated at the eastern, western faces, and $v$ is evaluated at the corners. The V-cell, where $v$ is evaluated at the centre, $\eta$ and $H$ are evaluated at the southern, northern faces, and $u$ is evaluated at the corners. The distance between two consecutive cells of the same type is $\Delta x$, zonally, and $\Delta y$, meridionally. The indices $i$ and $j$ correspond to the $i$-th zonal cell and the $j$-th meridional cell counted in the South-North direction and in the West-East direction respectively.](figs/arakawaCgrid.svg)

The mesh in use is an Arakawa Staggered regular C-grid (Arakawa, 1966) as illustrated in Fig. 2. It
is composed of three distinct cells: the U-cell, the V-cell, and
the T-cell, where at the centres are the $u$, the
$v$ and the $\eta$ variables of equations above.

The C-grid provides better precision for the non-linear advecting
terms than the B-grid, however it loses precision when evaluating
the Coriolis term in equations (Arakawa, 1966). For simplicity, the mesh will have constant
step-sizes $\Delta x$ and $\Delta y$. The indices $i$ and $j$ as
shown in Fig. 2 and in
Fig. 3 correspond to the $i$-th zonal cell
and the $j$-th meridional cell counted in the South-North
direction and in the West-East direction respectively.

### Boundary Conditions

Currently, only Dirichelet conditions are implemented at the
boundaries. Indeed, if the T-cells domain has $M\times N$ nodes
then the U-cells have $M\times (N+1)$ nodes and the V-cells have
$(M+1)\times N$ nodes. $\eta$ is calculated within
$\{2,\,...,\,(M-1)\} \times \{2,\,...,\,(N-1)\}$ and $u$ and $v$
are calculated within $\{2,\,...,\,(M-1)\}\times \{2,\,...,\,N\}$
and $\{2,\,...,\,M\}\times \{2,\,...,\,(N-1)\}$.

#### Null-flux

A land mask, $m_T$, for the T-cells mesh is introduced. The goal is to impose a null-flux boundary condition surrounding any land cell, i.e. 

$$
\vec{v} \cdot \vec{n} = 0.
$$

It returns 1 if the cell is filled with water and 0 if the cell is land. This implies the definition of appropriate null-fluxes masks, $m_U$ and $m_V$, for the U and V-cells. Thus, for every $i,j$ such that $m_T = 0$, it is required that $m_U = 0$, $m_{U\,i+1} = 0$, $m_V = 0$ and $m_{V\,j+1} = 0$. Everywhere else the value of the masks is $1$. The T, U and V masks are to be applied in the numerical scheme to the T-cell properties, the U-cell properties and the V-cell properties, respectively.

#### No-slip

The optional no-slip boundary condition (Pedlosky, 1987) consists of both null-flux and null-tangential velocities at the vertical walls of the domain, i.e.

$$
\vec{v} \cdot \vec{n} = 0,
$$

and

$$
\vec{v} \; \bot \; \vec{n} = 0.
$$

Thus, for every $i,j$ such that $m_T = 0$, it is required, additionally to the defined above null-flux condition, that $m_{U\,i,\,j+1} = 0$, $m_{U\,i+1,\,j+1} = 0$, $m_{U\,i,\,j-1} = 0$ and $m_{U\,i+1,\,j-1} = 0$ and that $m_{V\,i-1,\,j} = 0$, $m_{V\,i+1,\,j} = 0$, $m_{V\,i-1,\,j+1} = 0$ and $m_{V\,i+1,\,j+1} = 0$.

One interesting aspect of the no-slip boundary condition is that it necessarily requires a global zero-curl for closed domains,

$$
\oint \vec{v} \cdot \vec{dS} = 0.
$$

Hence, using the Kelvin-Stokes theorem, the no-slip boundary condition is an interesting configuration to test the correct implementation of the model: the curl within the domain must sum up to zero. Nevertheless, the no-slip boundary condition is a very strong constraint that acts on the kinematics and not on the dynamics of the motion per se (it is independent of the equation of motion).

**Fig. 3:** ![Detailed mesh emphasizing the boundaries. Composite of T, U and V-cells, the mesh illustrates the zone of integration of each type of cell: the blue rectangle contains the T-cells computed nodes, the thin green rectangle contains the U-cells computed nodes, the thin red rectangle contains the V-cells computed nodes. The thick green and red rectangles, however, delimit respectively the faces of the U and V-cells computed nodes.](figs/Boundaries2.svg)

#### Radiative Boundary Conditions

When no wall is to be considered at the boundaries, then all perturbation generated inside the domain eventually needs to go out of the domain. Furthermore, it could be interesting to propagate perturbations and information coming from outside of the domain. To this purpose are considered the broad class of open boundary conditions (OBC). The OBC are classified into two functional groups: the passive boundary conditions and the active boundary conditions. The passive boundary condition are designed to let information generated inside the domain to leave the domain, whereas the active boundary conditions try to propagate information from outside into the domain. 

Most regional oceanic modellers desire both aspects, of letting information out of and into the domain, which is considered a challenge. For very good reviews on OBC for regional ocean models, refer to Blayo and Debreu (2005) and Herzfeld (2008). For more recent types of radiative boundary conditions suitable for internal waves as well as the external mode, see Marsaleix et al. (2009). 

Radiative boundary conditions are passive boundary conditions (designed to let perturbations go out of the boundary) and usually consider the linearized hyperbolic version of the equations along the normal axis relatively to the open boundary. In this work, the gravity wave radiative method (also known as Sommerfeld radiative method) was implemented for the water elevation, $\eta$, and for the velocity tangential to the open boundary:

$$
\frac{\partial \Phi}{\partial t} + \vec{c} \cdot \vec{n} \frac{\partial \Phi}{\partial \vec{n}} = 0,
$$

where $\Phi$ is either the water elevation or the tangential velocity, $\vec{n}$ is the external normal vector to the open boundary and $\vec{c}$ is the phase wave celerity vector. In every occurence, the normal celerity wave intensity is considered to be $\vec{c} \cdot \vec{n} = \sqrt{g\,H}$. The passive Flather (1976) radiation method was implemented for the velocity normal to the open boundary:

$$
H\,\vec{v} \cdot \vec{n} = \eta \, \vec{c} \cdot \vec{n},
$$

where $\vec{v}$ is the flow velocity vector. Both methods are implemented with the normal velocity outside of the elevation node (NVOE).

The NVOE indicates how the radiative condition is implemented in the numerical scheme and it is dependent on the design of the grid and of the boundaries by the modeller. The configuration with the normal velocity inside the elevation node (NVIE) is not considered in this work, though it is thought to be relatively easier to adapt it from a NVOE configuration, rather than the other way around (Herzfeld, 2008). 

Finally, the NVIE was seen to return poorer results in some basic experiments (Herzfeld, 2008). The NVIE-NVOE dichotomy is pertinent as each implementation will affect differently each term of the SWE equations. Herzfeld (2008) reported some very interesting tables describing which terms of the SWE equations are affected by the NVIE and NVOE implementations.

### The Numerical Scheme

For simplicity in the notation, the indices $i$ and $j$ will be
omitted by default. The spatial finite difference numerical scheme
is centered in time and centered in time (CTCS) described in Kantha and Clayson (2000):

For the zonal momentum (U-Cell), the first-order spatial discretization writes:

$$
\begin{aligned}
\frac{\partial Hu}{\partial t} &= - \frac{\left(H u u\right)_{i + 1/2} - \left(H u u\right)_{i - 1/2}}{\Delta x} \\
& - \frac{\left(H u v\right)_{j + 1/2} - \left(H u v\right)_{j - 1/2}}{\Delta y} \\
& + f \left(H v\right) \\
& + \frac{\nu}{\Delta x} \left( \left(H\frac{\partial u}{\partial x}\right)_{i+1/2} - \left(H\frac{\partial u}{\partial x}\right)_{i-1/2} \right) \\
& + \frac{\nu}{\Delta y} \left( \left(H\frac{\partial u}{\partial y}\right)_{j+1/2} - \left(H\frac{\partial u}{\partial y}\right)_{j-1/2} \right) \\
& - \frac{gH(\eta _{i + 1/2} - \eta _{i - 1/2})}{\Delta x} \\
& + \frac{\rho _a}{\rho _0}C_a u_{10} \sqrt{u_{10}^2 + v_{10}^2} \\
& - C_D u_b \sqrt{u_b^2 + v_b^2} \\
& \equiv Ru
\end{aligned}
$$

where the halved indices correspond to fluxes at the U-cells'
faces. Thus, the CTCS fluxes write:

$$
\begin{aligned}
\left(H u u\right)_{i + 1/2} &= m_{U\,i+1}\,H(u_{i + 1} + u)^2 /4 \\
\left(H u u\right)_{i - 1/2} &= m_{U\,i-1}\,H_{i - 1} (u + u_{i - 1})^2 /4 \\
\left(H u v\right)_{j + 1/2} &= m_{U\,j+1}\,(H_{i - 1} + H_i + H_{i - 1,j + 1} + H_{i,j + 1}) \\
&\quad \times (u_{j + 1} + u)(v_{i - 1,j + 1} + v_{i,j + 1})/16 \\
\left(H u v\right)_{j - 1/2} &= m_{U\,j-1}\,(H_{i - 1} + H_i + H_{i - 1,j - 1} + H_{i,j - 1}) \\
&\quad \times (u + u_{j - 1})(v_{i - 1} + v)/16 \\
f\left(H v\right) &= f(H + H_{i - 1})/2 \\
&\quad \times \frac{(m_V\,v)_{i - 1} + m_V\,v + (m_V\,v)_{j + 1} + (m_V\,v)_{i - 1,j + 1}}{ m_{V\,i-1} + m_{V} + m_{V\,j+1} + m_{V\,i-1,j+1}}
\end{aligned}
$$

$$
\begin{aligned}
\nu \left(\left(H\frac{\partial u}{\partial x}\right)_{i+1/2} - \left(H\frac{\partial u}{\partial x}\right)_{i-1/2} \right) &= \nu \left( m_{U\,i+1} H \frac{u_{i+1} - u}{\Delta x} - m_{U\,i-1} H_{i-1} \frac{u-u_{i-1}}{\Delta x} \right) \\
\nu \left(\left(H\frac{\partial u}{\partial y}\right)_{j+1/2} - \left(H\frac{\partial u}{\partial y}\right)_{j-1/2} \right) &= \nu \left( m_{U\,j+1} H \frac{u_{j+1} - u}{\Delta y} - m_{U\,j-1} H_{j-1} \frac{u-u_{j-1}}{\Delta y} \right) \\
gH(\eta _{i + 1/2} - \eta _{i - 1/2}) &= g(H + H_{i - 1})/2 (\eta - \eta _{i - 1}) \\
C_D u\sqrt{u^2 + v^2} &= C_D u\sqrt{u^2 + \left(\frac{(m_V\,v)_{i - 1} + m_V\,v + (m_V\,v)_{j + 1} + (m_V\,v)_{i - 1,j + 1}}{ m_{V\,i-1} + m_{V} + m_{V\,j+1} + m_{V\,i-1,j+1}}\right)^2}
\end{aligned}
$$

Notice how the $\left(H u u\right)_{j + 1/2}$, $\left(H u v\right)_{j - 1/2}$, $f\left(H v\right)$ and $C_D u\sqrt{u^2 + v^2}$ terms lose significant precision compared to the other terms, due to their 4-term averaging.

Hence, rewriting the full momentum CTCS spatial scheme we get:

$$
\begin{align}
\frac{\partial Hu}{\partial t} = & - \frac{1}{\Delta x}\left( m_{U\,i+1} H(u_{i + 1} + u)^2 /4 - m_{U\,i-1} H_{i - 1} (u + u_{i - 1})^2 /4 \right) \\
& - \frac{1}{\Delta y}\biggl( m_{U\,j+1}(H_{i - 1} + H_i + H_{i - 1,j + 1} + H_{i,j + 1}) \\
& \quad \times (u_{j + 1} + u)(v_{i - 1,j + 1} + v_{i,j + 1})/16 \\
& \quad - m_{U\,j-1}(H_{i - 1} + H_i + H_{i - 1,j - 1} + H_{i,j - 1}) \\
& \quad \times(u + u_{j - 1})(v_{i - 1} + v)/16 \biggr) \\
& + f(H + H_{i - 1})/2 \\
& \quad \times \frac{(m_V\,v)_{i - 1} + m_V\,v + (m_V\,v)_{j + 1} + (m_V\,v)_{i - 1,j + 1}}{m_{V\,i-1} + m_{V} + m_{V\,j+1} + m_{V\,i-1,j+1}} \\
& + \frac{\nu}{\Delta x} \left( m_{U\,i+1} H \frac{u_{i+1} - u}{\Delta x} - m_{U\,i-1} H_{i-1} \frac{u-u_{i-1}}{\Delta x} \right) \\
& + \frac{\nu}{\Delta y} \left( m_{U\,j+1} H \frac{u_{j+1} - u}{\Delta y} - m_{U\,j-1} H_{j-1} \frac{u-u_{j-1}}{\Delta y} \right) \\
& - \frac{g(H + H_{i - 1})}{2\Delta x}(\eta - \eta _{i - 1}) \\
& + \frac{\rho _a}{\rho _0}C_a u_{10} \sqrt{u_{10}^2 + v_{10}^2} \\
& - C_D u\sqrt{u^2 + \left( \frac{(m_V\,v)_{i - 1} + m_V\,v + (m_V\,v)_{j + 1} + (m_V\,v)_{i - 1,j + 1}}{m_{V\,i-1} + m_{V} + m_{V\,j+1} + m_{V\,i-1,j+1}} \right)^2} \\
& \equiv Ru
\end{align}
$$

For the meridional spatial momentum scheme in the V-Cells, clever
symmetry one-to-one relations with zonal momentum scheme in the
U-cells are used:

- switch $\Delta x$ and $\Delta y$: $\Delta x \leftrightarrow \Delta y$
- switch $i$ and $j$: $i \leftrightarrow j$
- switch $u$ and $v$: $u \leftrightarrow v$
- switch signal of the Coriolis term: $\left( + \leftrightarrow - \right)$
- switch $M$ and $N$: $M \leftrightarrow N$

The finite-difference first-order numerical scheme for the waterlevel (T-Cell) writes out:

$$
\begin{aligned}
\frac{\partial \eta}{\partial t} &= - \frac{\left(Hu\right)_{i + 1/2} - \left(Hu\right)_{i - 1/2}}{\Delta x} \\
& - \frac{\left(Hv\right)_{j + 1/2} - \left(Hv\right)_{j - 1/2}}{\Delta y} \\
& \equiv R\eta
\end{aligned}
$$

and each face's CTCS flux term writes down:

$$
\begin{aligned}
\left(Hu\right)_{i + 1/2} &= m_{T\,i+1}\,(H + H_{i + 1})/2\;u_{i + 1} \\
\left(Hu\right)_{i - 1/2} &= m_{T\,i-1}\,(H_{i - 1} + H)/2\;u \\
\left(Hv\right)_{j + 1/2} &= m_{T\,j+1}\,(H + H_{j + 1})/2\;v_{j + 1} \\
\left(Hv\right)_{j - 1/2} &= m_{T\,j-1}\,(H_{j - 1} + H)/2\;v
\end{aligned}
$$

Thus, the full waterlevel CTCS numerical scheme is:

$$
\begin{aligned}
\frac{\partial \eta}{\partial t} &= - \frac{m_{T\,i+1}\,(H + H_{i + 1})/2\;u_{i + 1} - m_{T\,i-1}\,(H_{i - 1} + H)/2\;u}{\Delta x} \\
& - \frac{m_{T\,j+1}\,(H + H_{j + 1})/2\;v_{j + 1} - m_{T\,j-1}\,(H_{j - 1} + H)/2\;v}{\Delta y} \\
& \equiv R\eta
\end{aligned}
$$

The time scheme used is the Leapfrog as described in Kantha and Clayson (2000):

$$
\begin{aligned}
\eta ^{l + 1} &= \eta ^{l - 1} + 2\Delta t\,R\eta \\
H^{l + 1} &= \eta ^{l + 1} + d \\
u^{l + 1} &= \frac{u^{l - 1} (H^{l - 1} + H_{i - 1}^{l - 1}) + 4\Delta t\,Ru}{H^{l + 1} + H_{i - 1}^{l + 1}} \\
v^{l + 1} &= \frac{v^{l - 1} (H^{l - 1} + H_{j - 1}^{l - 1}) + 4\Delta t\,Rv}{H^{l + 1} + H_{j - 1}^{l + 1}}
\end{aligned}
$$

Notice how the leapfrog time scheme obliges two initial conditions
at $t_0$ and at $t_1$. Hence, in order to avoid mode decoupling, a
Robert-Asselin filter (Asselin, 1972) for $u,\,v,\,\eta$ at each integration time-step is used,
as suggested by Kantha and Clayson (2000):

$$P^l = P^l + \gamma \left(P^{l - 1} - 2P^l + P^{l + 1}\right)$$

where $\gamma$ is a parameter set to $0.1$ (Kantha and Clayson, 2000). The Robert-Asselin
provides a good coupling between the two initial conditions, at the expense of some loss in precision (Asselin, 1972).

The radiative scheme implemented follows a NVOE stencil on a C grid (Herzfeld, 2008). The *western* boundary radiative condition is defined, for the elevation and the component of velocity perpendicular to the boundary, by:

$$
\begin{aligned}
\eta^{l+1}_{1,\,j} &= \eta_{1,\,j} - 2 \frac{\Delta t}{\Delta x} \,\sqrt{g\, H_{1,\,j}} 
\, \left(\eta_{1,\,j} - \eta_{2,\,j}\right) \\
u^{l+1}_{1,\,j} &= - \sqrt{\frac{g}{H^{l+1}_{1,\,j}}} \, \eta^{l+1}_{1,\,j}
\end{aligned}
$$

for $j=1,\,...,\,N$, and is defined by, for the velocity component tangent to the boundary:

$$
v^{l+1}_{1,\,j}= \frac{v_{1,\,j} \, \left(H_{1,\,j} + H_{1,\,j-1}\right) - 2 \, \frac{\Delta t}{\Delta x} \, \sqrt{g\,\frac{H_{1,\,j} + H_{1,\,j-1}}{2}} \, \left(v_{1,\,j} - v_{2,\,j}\right)}{H^{l+1}_{1,\,j} + H^{l+1}_{1,\,j-1}}
$$

for $j=2,\,...,\,N$.

For the *eastern* boundary, the radiation boundary condition writes:

$$
\begin{aligned}
\eta^{l+1}_{M,\,j} &= \eta_{M,\,j} - 2 \frac{\Delta t}{\Delta x} \,\sqrt{g\, H_{M,\,j}} 
\, \left(\eta_{M,\,j} - \eta_{M-1,\,j}\right) \\
u^{l+1}_{M+1,\,j} &= - \sqrt{\frac{g}{H^{l+1}_{M,\,j}}} \, \eta^{l+1}_{M,\,j}
\end{aligned}
$$

for $j=1,\,...,\,N$, and is defined by, for the velocity component tangent to the boundary:

$$
v^{l+1}_{M,\,j}= \frac{v_{M,\,j} \, \left(H_{M,\,j} + H_{M,\,j-1}\right) - 2 \, \frac{\Delta t}{\Delta x} \, \sqrt{g\,\frac{H_{M,\,j} + H_{M,\,j-1}}{2}} \, \left(v_{M,\,j} - v_{M-1,\,j}\right)}{H^{l+1}_{M,\,j} + H^{l+1}_{M,\,j-1}}
$$

for $j=2,\,...,\,N$. Note that the Flather (1976) radiation condition applied to the normal component of velocity to the boundary can be replaced with a simple null-gradient and yield similar results:

$$
u^{l+1}_{1,\,j} = u^{l+1}_{2,\,j}
$$

for $i=1$, and

$$
u^{l+1}_{m+1,\,j} = u^{l+1}_{m,\,j}
$$

for $i=m+1$.

Once more, to derive an adequate scheme for the *southern* and *northern* boundary conditions, simply follow the symmetrical rules below and apply them to the preceding equations:

- switch $i$ and $j$: $i \leftrightarrow j$
- switch $u$ and $v$: $u \leftrightarrow v$
- switch $M$ and $N$: $M \leftrightarrow N$
- switch $\Delta x$ and $\Delta y$: $\Delta x \leftrightarrow \Delta y$
- switch (*West*, *East*) with (*South*, *North*)

The stability criterion is the Courant-Friedrich-Levy criterion (Courant et al., 1959) described in
Kantha and Clayson (2000):

$$
\Delta t \left( \sqrt{gH} + V_{max} \right) \left( \frac{1}{\Delta x} + \frac{1}{\Delta y} \right) < 1.0
$$

where $V_{max}$ is the maximum advection field intensity in m/s.
Note that for stability reasons, in the momentum equations, the
friction terms are evaluated at time $l-1$.
