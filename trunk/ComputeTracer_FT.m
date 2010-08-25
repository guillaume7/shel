%--------------------------------------------------------------------------
%     This file is part of SHEL SHallow-water numerical modEL
% 
%     SHEL is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     SHEL is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with SHEL.  If not, see <http://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------
function [Tr, dttr] = ComputeTracer_FT(Tr, K_L)
%Advection-diffusion tracer equation solver
%
%Explicit forward in time and upwind
%
%         Arakawa C staggered grid (Arakawa 1966)
%
%    --j y
%   |           --- V ---           V ------- V         U -- T -- U
%   i          |         |          |         |         |         |
%   x          U    T    U          T    U    T         |    V    |
%              |         |          |         |         |         |
%               --- V ---           V ------- V         U -- T -- U
%                 T-cell               U-cell              V-cell
%
%       T(eta, H, tr, d, f, mask,   U(u, mask_u         V(v, mask_v
%           x, y, wind and bottom      x_u, y_u)           x_v, y_v)
%            stress)
%
global dx;
global dy;
global dt;
 
%T-cells
global H;
global H_new;
global mask;
 
%U-cells
global u_a;
global mask_u;

%V-cells
global v_a;
global mask_v;

%% Solve the tracer spatial scheme in the whole of the domain
    RHSt = mask .* ( ...
        ComputeSpaceT_UP_U(H,u_a,mask_u, Tr, K_L, dx) ... %Along U
       + ComputeSpaceT_UP_U(H',v_a',mask_v',Tr',K_L,dy)' ... %then along V
    );

%% FIND THE OPTIMAL TRACER DT AS A MULTIPLE OF MOMENTUM DT
% within the reasonable limit of maximum 100 * dt.
    minRHSt  = min(min(RHSt));
    ind = find(RHSt == minRHSt);
    dttr = min( floor( abs(- H(ind) .* Tr(ind) ./ minRHSt) / dt ), 1000) ...
        * dt;

%% TRACER FORWARD-IN-TIME SCHEME
    Tr_new = mask .* (H .* Tr + dttr * RHSt) ./ H_new;

%% Update matrices
    Tr = Tr_new;

function RHSt = ComputeSpaceT_UP_U(H, u, mask_u, Tr, K_L, dx_L)
%function RHSt = ComputeSpaceT_UP_U(H, u, mask_u, Tr, K, dx)
%
%UPWIND
%
    global nullgrad;
    nullgrad = false;
    
    [M,N]=size(H);
    Hm = zeros(M+1,N);
    TrUm = zeros(M+1,N);
    TrUp = zeros(M+1,N);
    TrDif = zeros(M+1,N);

    %Computing quantities for the interior of the domain
    TrUm(1:M,:) = .5 * (abs(u(1:M,:)) - u(1:M,:)) .* Tr;
    TrUp(2:M+1,:) = .5 * (abs(u(2:M+1,:)) + u(2:M+1,:)) .* Tr;
    TrDif(2:M,:) = Tr(1:M-1,:) - Tr(2:M,:);
    Hm(2:M,:) = .5 * ( H(1:M-1,:) + H(2:M,:) );
    
    %Computing the quantities for the boundaries
    if nullgrad %Considering null-gradient:
        TrUm(M+1,:) = .5 * (abs(u(M+1,:)) - u(M+1,:)) .* Tr(M,:);
        TrUp(1,:) = .5 * (abs(u(1,:)) + u(1,:)) .* Tr(1,:);
        TrDif(1:M:M+1,:) = 0.;
    else %else null value is considered
        TrUm(M+1,:) = .5 * (abs(u(M+1,:)) - u(M+1,:)) * 0.;
        TrUp(1,:) = .5 * (abs(u(1,:)) + u(1,:)) * 0.;
        TrDif(1,:) = - Tr(1,:);
        TrDif(M+1,:) = Tr(M,:);
    end
    
    %This algo considers that the bathymetry has always a null gradient
    %at the boundary:
    Hm(1:M:M+1,:) = H(1:M-1:M,:);

    %Computing the fluxes at each face of the M+1 faces
    FFluxU = ComputeFaceFluxT_UP(Hm, TrUp, TrUm, TrDif, dx_L, K_L);
    
    %The T-cell increment is the flux balance.
    RHSt = mask_u(1:M,:) .* FFluxU(1:M,:) ...
         - mask_u(2:M+1,:) .* FFluxU(2:M+1,:);

%Pure upwind scheme.
function RHSt = ComputeFaceFluxT_UP(Hm, TrUp, TrUm, TrDif, dx_L, K_L)
    RHSt = Hm .* ( ...
            TrUp - TrUm ...
            + K_L / dx_L * TrDif ...
    ) / dx_L;

%--------------------------------------------------------------------------
%SHEL SHallow-water numerical modEL
%Copyright (C) 2006,2009,2010. Guillaume Riflet, Instituto Superior Técnico
%da Universidade Técnica de Lisboa.
%--------------------------------------------------------------------------
