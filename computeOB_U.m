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
function [eta_new, H_new, u_new, v_new] = computeOB_U( ...
    eta_new, eta_old, H_new, H_old, d,...
    u_new, mask_u, ...
    v_new, v_old, mask_v, ...
    dt, dx ...
)
%function computeOB
%Solve the north, south, east and west boundary condition
%for the eta, u and v variables.

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

global g; 
  
%BC conditions 
global flather;
global neumann;
global radiatetan; 
global radiatelevel;

flather = 0;
neumann = 0;
radiatetan = 0; 
radiatelevel = 1;

[M,N] = size(d);

%% Eastern (1,:) and Western (M,:) boundaries    
if radiatelevel 
    %Waterlevel: GW radiation 
    eta_new(1:M-1:M,:) = mask_u(1:M:M+1,:) .* ( ...
        eta_old(1:M-1:M,:)- dt / dx .* sqrt( g * H_old(1:M-1:M,:) ) ... 
        .* ( eta_old(1:M-1:M,:) - eta_old(2:M-3:M-1,:) ) ...
    ); 
    %no mask below, cuz H must yield -99 where land is.W
    H_new(1:M-1:M,:) = eta_new(1:M-1:M,:) + d(1:M-1:M,:);
end 
if flather 
    %Normal velocity: Flather 
    u_new(1:M:M+1,:) = [-1 * ones(1,N); ones(1,N)] .*  mask_u(1:M:M+1,:) ...
        .* sqrt( g ./ H_new(1:M-1:M,:) ) .* eta_new(1:M-1:M,:); 
elseif neumann
    %Normal velocity: Null-gradient 
    u_new(1:M:M+1,1:N) = mask_u(1:M:M+1,1:N) .* u_new(2:M-2:M,1:N); 
end 
if radiatetan 
    %Tangential velocity: GW radiation 
    v_new(1:M-1:M,2:N) = mask_v(1:M-1:M,2:N) .* ... 
        ( ... 
            v_old(1:M-1:M,2:N) .* ( H_old(1:M-1:M,2:N) + H_old(1:M-1:M,1:N-1) ) ... 
            - dt / dx .* sqrt( g * .5 * ( H_old(1:M-1:M,2:N) + H_old(1:M-1:M,1:N-1) ) ) ... 
            .* ( v_old(1:M-1:M,2:N) - v_old( 2:M-3:M-1,2:N) ) ... 
        ) ... 
        ./ ( H_new(1:M-1:M,2:N) + H_new(1:M-1:M,1:N-1) ); 
end 

%--------------------------------------------------------------------------
%SHEL SHallow-water numerical modEL
%Copyright (C) 2006,2009,2010. Guillaume Riflet, Instituto Superior Técnico
%da Universidade Técnica de Lisboa.
%--------------------------------------------------------------------------
