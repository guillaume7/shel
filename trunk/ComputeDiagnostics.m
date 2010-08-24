%--------------------------------------------------------------------------
%     This file is part of SHEL SHallow-water numerical modEL
% 
%     SHEL is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     Foobar is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------
function diag = ComputeDiagnostics(l)
%function ComputeDiagnostics
%
%Computes global:
% volume (m3)
% kinetic energy (J)
% potential energy (J)
% total energy (J)
% vorticity (m2/s)
% enstrophy (J)
% squared normal strain (J)
% squared shear strain (J)
% squared strain (J)
% Okubo-Weiss
% momentum (m^3 * m / s)a
%
%Computes local:
% kinetic energy (J)
% potential energy (J)
% total energy (J)
% curl (1/s)
% divergence (1/s)
% stretch (normal strain) (1/s)
% shear (shear strain) (1/s)
% squared strain (1/s)
% enstrophy (1/s2)
% OkuboWeiss (1/s2)
%
%
%    --j y
%   |           --- V ---           V ------- V         U -- T -- U
%   i          |         |          |         |         |         |
%   x          U    T    U          T    U    T         |    V    |
%              |         |          |         |         |         |
%               --- V ---           V ------- V         U -- T -- U
%                 T-cell               U-cell              V-cell
%                  M x N               M+1 x N             M x N+1
%       T(eta, H, tr, d, f, mask,   U(u, mask_u         V(v, mask_v
%           x, y, wind and bottom      x_u, y_u)           x_v, y_v)
%            stress, curl(v,u,0))
%
%    --j y
%   |           --- U ---  
%   i          |         | 
%   x          V    W    V 
%              |         | 
%               --- U ---  
%                 T-cell
%                M+1 x N+1
%       W( curl(u,v,0) )
%        (deformation rate = curl(-u,v,0)
%%%

global u_a;
global v_a;
global u;
global v;
global u_t;
global v_t;
global H;
global eta;
global M;
global N;
global dx;
global dy;
global dt;
global rho0;
global g;
global K;
global mask;
global ke;
global pe;
global strechrate_t;
global shearrate_w;
global divergence_t;
global sqstrain_t;
global sqstrechrate_t;
global sqshearrate_w;
global curl_w;
global enstrophy_t;
global enstrophy_w;
global okuboweiss_t;

%Eddy kinetic energy
global gradux;
global graduy;
global gradvx;
global gradvy;
global eke;

%time-dependent global properties
global time;
global volume;
global iKe;
global iPe;
global vtime;
global iVort;
global iEnst;
global iSqStrech;
global iSqShear;
global iSqStrain;
global iOWeiss;
global iMomentumU;
global iMomentumV;
global iEke;

dA = dx * dy; %m2

%Computes the average velocity field
u_a = ((l - 1) * u_a + u ) / l;
v_a = ((l - 1) * v_a + v ) / l;

%Eddy kinetic energy
gradux = (u(2:M+1,:) - u(1:M,:)) / dx;
gradvy = (v(:,2:N+1) - v(:,1:N)) / dy;
graduy(2:M,2:N-2) = (u(2:M,3:N-1) - u(2:M,2:N-2)) / dy;
gradvx(2:M-2,2:N) = (v(3:M-1,2:N) - v(2:M-2,2:N)) / dx;

%Compute diagnostic fields
%Curl
computeCurl_w; % 1 / s

%Strain rate
strechrate_t = mask .* ((u(2:M+1,:) - u(1:M,:)) * dy ...
            - (v(:,2:N+1) - v(:,1:N)) * dx) / dA; % 1 / s
computeShearRate_w % 1 / s

%Growth rate
divergence_t = mask .* ((u(2:M+1,1:N) - u(1:M,1:N)) * dy ...
                + (v(1:M,2:N+1) - v(1:M,1:N)) * dx) / dA; % 1/ s

%Potential vorticity
computePotentialVorticity_w

%This is the so-called scalar of Arakawa-Okubo-Weiss.
%Computes the local enstrophy, squared strain rate and Okubo-Weiss
%parameter %Look for Weiss
%La Jolla Institute Report from 1981.
%Physica D: Nonlinear 1991.
%Look for Okubo-Weiss parameter/function
%Look for Arakawa 1966 paper.
enstrophy_w = .5 * ( curl_w.^2); %1 / s2
enstrophy_t = interpolFtoT(enstrophy_w, mask, M, N);

sqshearrate_w = .5 * shearrate_w.^2;
sqshearrate_t = interpolFtoT(sqshearrate_w, mask, M, N);

sqstrechrate_t = .5 * strechrate_t.^2 .* mask;

sqstrain_t = sqstrechrate_t + sqshearrate_t;

sqdivergence_t = .5 * divergence_t.^2 .* mask;

okuboweiss_w = (sqshearrate_w - enstrophy_w);
okuboweiss_t = interpolFtoT(okuboweiss_w, mask, M, N);
okuboweiss_t = okuboweiss_t + sqstrechrate_t .* mask;
%versão do Aires do escalar de Okubo-Weiss
okuboweiss_t = okuboweiss_t  - sqdivergence_t;

%%interpolate from U and V to T cells
u_t = .5 * ( u(1:M,:) + u(2:M+1,:) );
v_t = .5 * ( v(:,1:N) + v(:,2:N+1) );

Hnoland = H;
Hnoland(find(H < -1)) = 0.;

%Volume
%volume(l) = dA * sum(sum(Hnoland)); %m3
%Perturbed volume
volume(l) = dA * sum(sum(eta .* mask)); %m3

%Momentum (actually it's only the integrated velocity)
iMomentumU(l) = dA * sum(sum(u_t .* Hnoland,2)); % kg * m / s
iMomentumV(l) = dA * sum(sum(v_t .* Hnoland,1)); % kg * m / s

%Kinetic energy
ke = .5 * rho0 * dA * (u_t.^2 + v_t.^2) .* Hnoland .* mask; % kg * m * m / s2 = J
iKe(l) = sum(sum(ke)); % kg * m * m / s2 = J

%Potential energy
pe = .5 * rho0 * g * dA .* eta.^2 .* mask; %kg m * m / s2 = J
iPe(l) = sum(sum(pe)); %kg m * m / s2 = J

%Eddy kinetic energy
eke = eke + rho0 * dA * K * 2 * dt * (gradux.^2 + gradvy.^2 + graduy.^2 + gradvx.^2) ...
    .* Hnoland .* mask;
iEke(l) = sum(sum(eke));

%global vorticity
iVort(l) = sum(sum(curl_w)) * dA; % m2 / s

%global enstrophy
iEnst(l) = sum(sum(enstrophy_w)) * dA; % m2 / s2

%global squared strain
iSqShear(l) = sum(sum(sqshearrate_w))* dA;
iSqStrech(l) = sum(sum(sqstrechrate_t)) * dA;
iSqStrain(l) = sum(sum(sqstrain_t)) * dA;

%global Okubo-Weiss
iOWeiss(l) = sum(sum(okuboweiss_t)) * dA; %<--- too large numerical error
%iOWeiss(l) = iSqStrain(l) - iEnst(l);

%time
vtime(l) = time;

diag = [...
        volume(l), ...
        iKe(l), ...
        iPe(l), ...
        iVort(l), ...
        iEnst(l), ...
        iSqStrain(l), ...
        iOWeiss(l) ...
        ];

%% ComputeCurl_w % 1 / s
function computeCurl_w
%function computeCurl_w

global curl_w;
global curl_t;
global mask;
global mask_v;
global mask_u,
global u;
global v;
global M;
global N;
global dy;
global dx;

dA = dx .* dy;

%Curl of velocity in the interior...
curl_w(2:M,2:N) = ((mask_v(2:M,2:N) .* v(2:M,2:N) ...
    - mask_v(1:M-1,2:N) .* v(1:M-1,2:N)) .* dy ...
    - (mask_u(2:M,2:N) .* u(2:M,2:N) ...
    - mask_u(2:M,1:N-1) .* u(2:M,1:N-1)) .* dx) ./ dA; % 1 / s
%... at the corners %Thank you for the conversation with F.L.!
%lower-left SW
curl_w(1,1) = ((mask_v(1,1) * v(1,1)) * dy ...
    - (mask_u(1,1) * u(1,1)) * dx)/dA; % 1 / s
%lower-right SE
curl_w(1,N+1) = ((mask_v(1,N+1) * v(1,N+1)) * dy ...
    - ( - mask_u(1,N) * u(1,N)) * dx)/dA; % 1 / s
%top-left NW
curl_w(M+1,1) = (( - mask_v(M,1) * v(M,1)) * dy ...
    - (mask_u(M+1,1) * u(M+1,1)) * dx)/dA; % 1 / s
%top-right NE
curl_w(M+1,N+1) = ((- mask_v(M,N+1) * v(M,N+1)) * dy ...
    - (- mask_u(M+1,N) * u(M+1,N)) * dx)/dA; % 1 / s
%... at the western boundary
curl_w(2:M,1) = ((mask_v(2:M,1) .* v(2:M,1) - mask_v(1:M-1,1) .* v(1:M-1,1)) * dy ...
    - (mask_u(2:M,1) .* u(2:M,1)) * dx)/dA; % 1 / s
%... at the eastern boundary
curl_w(2:M,N+1) = ((mask_v(2:M,N+1) .* v(2:M,N+1) - mask_v(1:M-1,N+1) .* v(1:M-1,N+1)) * dy ...
    - (- mask_u(2:M,N) .* u(2:M,N)) * dx)/dA; % 1 / s
%... at the southern boundary
curl_w(1,2:N) = ((mask_v(1,2:N) .* v(1,2:N)) * dy ...
    - (mask_u(1,2:N) .* u(1,2:N) - mask_u(1,1:N-1) .* u(1,1:N-1)) * dx)/dA; % 1 / s
%... at the northern boundary
curl_w(M+1,2:N) = ((- mask_v(M,2:N) .* v(M,2:N)) * dy ...
    - (mask_u(M+1,2:N) .* u(M+1,2:N) - mask_u(M+1,1:N-1) .* u(M+1,1:N-1)) * dx)/dA; % 1 / s
% interpolation from the {Z,W,F}-cell to the T-Cell
curl_t = interpolFtoT(curl_w, mask, M, N);

%% ComputeShearRate_w % 1 / s
function computeShearRate_w
%function computeShearRate_w 1 / s
%

global shearrate_w;
global shearrate_t;
global mask;
global u;
global v;
global M;
global N;
global dy;
global dx;

dA = dx * dy;

%Curl of (-u,v,0)
%interior of the domain...
shearrate_w(2:M,2:N) = ((v(2:M,2:N) - v(1:M-1,2:N)) * dy ...
    + (u(2:M,2:N) - u(2:M,1:N-1)) * dx) / dA;
%bottom-left corner (SW)
shearrate_w(1,1) = ((v(1,1)) * dy ...
    + (u(1,1)) * dx) / dA;
%bottom-right corner (SE)
shearrate_w(1,N+1) = ((v(1,N+1)) * dy ...
    + (- u(1,N)) * dx) / dA;
%top-left corner (NW)
shearrate_w(M+1,1) = ((- v(M,1)) * dy ...
    + (u(M+1,1)) * dx) / dA;
%top-right corner (NE)
shearrate_w(M+1,N+1) = ((- v(M,N+1)) * dy ...
    + (- u(M+1,N)) * dx) / dA;
%western boundary
shearrate_w(2:M,1) = ((v(2:M,1) - v(1:M-1,1)) * dy ...
    + (u(2:M,1)) * dx) / dA;
%eastern boundary
shearrate_w(2:M,N+1) = ((v(2:M,N+1) - v(1:M-1,N)) * dy ...
    + (- u(2:M,N)) * dx) / dA;
%southern boundary
shearrate_w(1,2:N) = ((v(1,2:N)) * dy ...
    + (u(1,2:N) - u(1,1:N-1)) * dx) / dA;
%northern boundary
shearrate_w(M+1,2:N) = ((- v(M,2:N)) * dy ...
    + (u(M+1,2:N) - u(M+1,1:N-1)) * dx) / dA;
% interpolate from F-cell to a T-cell
shearrate_t = interpolFtoT(shearrate_w, mask, M, N);

%% ComputePotentialVorticity_w % 1 / s / m
function computePotentialVorticity_w
%function computePotentialVorticity_w
%

global potvorticity_t;
global mask;
global curl_w;
global potvorticity_w;
global coriolis;
global f;
global H;
global M;
global N;

%Adding coriolis
if coriolis
  potvorticity_w = curl_w + f; % ??? curl_w * dA  + f???
end

%potential vorticity
%within the domain...
potvorticity_w(2:M,2:N) = potvorticity_w(2:M,2:N) * 4. ./ ...
    (H(2:M,2:N) + H(1:M-1,2:N) + H(1:M-1,1:N-1) + H(2:M,1:N-1));
%bottom-left corner (SW)
potvorticity_w(1,1) = potvorticity_w(1,1) / ...
    (H(1,1));
%bottom-right corner (SE)
potvorticity_w(1,N+1) = potvorticity_w(1,N+1)/ ...
    (H(1,N));
%top-left corner (NW)
potvorticity_w(M+1,1) = potvorticity_w(M+1,1) / ...
    (H(M,1));
%top-right corner (NE)
potvorticity_w(M+1,N+1) = potvorticity_w(M+1,N+1) / ...
    (H(M,1));
%western boundary
potvorticity_w(2:M,1) = potvorticity_w(2:M,1) * 2. ./ ...
    (H(2:M,1) + H(1:M-1,1));
%eastern boundary
potvorticity_w(2:M,N+1) = potvorticity_w(2:M,N+1) * 2. ./ ...
    (H(1:M-1,N) + H(2:M,N));
%southern boundary
potvorticity_w(1,2:N) = potvorticity_w(1,2:N) * 2. ./ ...
    (H(1,2:N) + H(1,1:N-1));
%northern boundary
potvorticity_w(M+1,2:N) = potvorticity_w(M+1,2:N) * 2. ./ ...
    (H(M,2:N) + H(M,1:N-1));
%interpolate from an F-cell to a T-cell
potvorticity_t = interpolFtoT( potvorticity_w, mask, M, N);

function tgrid = interpolFtoT(fgrid, mask, M, N)
%% function tgrid = interpolFtoT(fgrid, mask, M, N)
tgrid = .25 * mask .* ( ...
              fgrid(1:M,1:N) ...
            + fgrid(2:M+1,1:N) ...
            + fgrid(2:M+1,2:N+1) ...
            + fgrid(1:M,2:N+1) ...
            );
%--------------------------------------------------------------------------
%SHEL SHallow-water numerical modEL
%Copyright (C) 2006,2009,2010. Guillaume Riflet, Instituto Superior Técnico
%da Universidade Técnica de Lisboa.
%--------------------------------------------------------------------------
