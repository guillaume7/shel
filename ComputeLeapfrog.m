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
function ComputeLeapfrog
%Shallow-water equations solver
%
%Explicit leapfrog scheme
%combined with Asselin-Roberts filter
%(the filter supresses the computational mode
%and prevents alternate-time-decoupling)
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
global gama;
global dx;
global dy;
global dt; 
 
%T-cells 
global eta_old;
global eta;
global eta_new; 
global H_old; 
global H;
global H_new; 
global d; 
global mask; 
 
%U-cells 
global u_old; 
global u;
global u_new; 
global mask_u; 
 
%V-cells 
global v_old; 
global v;
global v_new; 
global mask_v;

%% Solve the spatial schemes in the interior of the domain
    % For the water elevation
    RHSeta = ComputeContinuity(H, eta, u, v, mask);
    
    % For the u and v components of velocity
    % % The old 10%-faster way (but a lot more error-prone)
    %     [RHSu, RHSv] = ComputeUV2D_old;
    % Switch v <-> u; v_old <-> u_old; dy <-> dx and -1. <-> 1. (coriolis)
    % Transpose the matrices.
    RHSu = ComputeSpaceU_CS( H, H_old, eta, ...
                        u, u_old, v, v_old, ...
                        mask_u, dx, dy, 1.);
    RHSv = ComputeSpaceU_CS( H', H_old', eta', ...
                        v', v_old', u', u_old', ...
                        mask_v', dy, dx, -1.)';
 
%% LEAPFROG TIME-SCHEME
    % Water level time scheme
    eta_new = eta_old + 2 * dt * RHSeta .* mask; 
    H_new = eta_new + d;
    %u and v time schemes inside the domain. 
    %Leapfrog means that 2xdt is passed as teh time step.
    v_new = ComputeTimeU_FT(v_old', v_new', RHSv', H_old', H_new', mask_v', 2*dt)';
    u_new = ComputeTimeU_FT(u_old, u_new, RHSu, H_old, H_new, mask_u, 2*dt);

%% time scheme at the northern/southern and eastern/southern OB (2 x d
%% <- LeapFrog)

    %Southern and northern boundaries (along V -> transpose inputs)
    %(transposed inputs -> transposed outputs).
    [eta_new, H_new, v_new, u_new] = computeOB_U( ...
            eta_new', eta_old', H_new', H_old', d',...
            v_new', mask_v', ...
            u_new', u_old', mask_u', ...
            2*dt, dy ...
     );
    %Eastern and western boundaries (along U)
    %The previous outputs were transposed, so we transpose
    %them back as inputs.
    [eta_new, H_new, u_new, v_new] = computeOB_U( ...
            eta_new', eta_old, H_new', H_old, d,...
            u_new', mask_u, ...
            v_new', v_old, mask_v, ...
            2*dt, dx ...
     );
     
%% Asselin-Robert filter 
    eta = eta + gama * ( eta_old - 2 * eta + eta_new ); 
    u = u + gama * ( u_old - 2 * u + u_new ); 
    v = v + gama * ( v_old - 2 * v + v_new ); 

%% Update matrices
    eta_old = eta;
    eta = eta_new;
     
    H_old = H;
    H = H_new;
 
    u_old = u;
    u = u_new;
     
    v_old = v;
    v = v_new;

%--------------------------------------------------------------------------
%SHEL SHallow-water numerical modEL
%Copyright (C) 2006,2009,2010. Guillaume Riflet, Instituto Superior Técnico
%da Universidade Técnica de Lisboa.
%--------------------------------------------------------------------------
