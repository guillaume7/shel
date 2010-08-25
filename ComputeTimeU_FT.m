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
function u_new = ComputeTimeU_FT(u_old, u_new, RHSu, H_old, H_new, mask_u, dt)
%function u_new = ComputeTimeU_FT(u_old, u_new, RHSu, H_old, H_new, mask_u, dt)        
%
%Numerical solver with the Forward in Time (FT) Euler time stepping method.
%
    [M,N] = size(H_old);
    
    u_new(2:M,2:N-1) = mask_u(2:M,2:N-1) .* ( ... 
                u_old(2:M,2:N-1) .* ( H_old(2:M,2:N-1) + H_old(1:M-1,2:N-1) ) ... 
                + 2 * dt * RHSu(2:M,2:N-1) ... 
              ) ... 
              ./ ( H_new(2:M,2:N-1) + H_new(1:M-1,2:N-1) );    
          
%--------------------------------------------------------------------------
%SHEL SHallow-water numerical modEL
%Copyright (C) 2006,2009,2010. Guillaume Riflet, Instituto Superior Técnico
%da Universidade Técnica de Lisboa.
%--------------------------------------------------------------------------
