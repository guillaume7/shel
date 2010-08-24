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
function ExteriorReference
%function status = ExteriorReference
%
%Use this function when the boundaries are not blosed (i.e. they're open)
%to define the exterior reference solution.
%
%The numerical model uses the outer perimeter cells to read any imposed
%exterior solution (i=1,M && j=1,N).
%
%So far, available, are the classic null condition (Dirichelet, clamped)
% and the null-gradient condition (Neumann).

global Neumann;

Neumann = false;

%% Setup the reference condition at the open boundaries %%%%%%%%%%%%%%%%%%%%
if Neumann   
    
    Neumann_f
    
else
    
    Dirichelet_f
    
end

% End of setup of open boundaries %%%%%%%%%%%%%%%%%%%

function Dirichelet_f
%% Setup the Dirichelet (clamped) open boundaries %%%%%%%%%%%%%%%%%%%%

global d;
global H;
global H_old; %to use in dissipation terms
global eta;
global eta_old;
global u;
global u_old; %to use in dissipation terms
global v;
global v_old;
global N;
global M;

    
%% West i=1 %%
    
    %T-cells
    for j=2:N-1
        
        eta(1,j) = 0.;
        eta_old(1,j) = 0.;
        
        H(1,j) = d(1,j);    
        H_old(1,j) = d(1,j);    
        
    end

    %U-cells
    for j=2:N-1
        
        u(1,j) = 0.;
        u_old(1,j) = 0.;
        
    end
    
    %V-cells
    for j=2:N
        
        v(1,j) = 0.;
        v_old(1,j) = 0.;
        
    end

%% North j=1 %%
    
    %T-cells
    for i=2:M-1
        
        eta(i,1) = 0.;
        eta_old(i,1) = 0.;
        
        H(i,1) = d(i,1);    
        H_old(i,1) = d(i,1);    
        
    end
    
    %U-cells
    for i=2:M
        
        u(i,1) = 0.;
        u_old(i,1) = 0.;
        
    end

    %V-cells
    for i=2:M-1
        
        v(i,1) = 0.;
        v_old(i,1) = 0.;

    end
    
%% East i=M+1(U) || M(T,V) %%
    
    %T-cells
    for j=2:N-1
        
        eta(M,j) = 0.;
        eta_old(M,j) = 0.;
        
        H(M,j) = d(M,j);    
        H_old(M,j) = d(M,j);    
                
    end
    
    %U-cells
    for j=2:N-1
        
        u(M+1,j) = 0.;
        u_old(M+1,j) = 0.;
        
    end
    
    %V-cells
    for j=2:N
                
        v(M,j) = 0.;
        v_old(M,j) = 0.;

    end

%% South j=N(T,U) OR N+1(V) %%

    %T-cells
    for i=2:M-1

        eta(i,N) =0.;
        eta_old(i,N) = 0.;
        
        H(i,N) = d(i,N);    
        H_old(i,N) = d(i,N);    
        
    end

    %U-cells
    for i=2:M
                
        u(i,N) = 0.;
        u_old(i,N) = 0.;
        
    end
    
    %V-cells
    for i=2:M-1
                        
        v(i,N+1) = 0.;
        v_old(i,N+1) = 0.;

    end
% End of setup of Dirichelet open boundaries %%%%%%%%%%%%%%%%%%%

function Neumann_f
%% Setup the Neumann open boundaries %%%%%%%%%%%%%%%%%%%%
    
global H;
global H_old; %to use in dissipation terms
global eta;
global eta_old;
global u;
global u_old; %to use in dissipation terms
global v;
global v_old;
global N;
global M;

    
%% West i=1 %%
    
    %T-cells
    for j=2:N-1
        
        eta(1,j) = eta(2,j);
        eta_old(1,j) = eta_old(2,j);
        
        H(1,j) = H(2,j);    
        H_old(1,j) = H_old(2,j);    
        
    end

    %U-cells
    for j=2:N-1
        
        u(1,j) = u(2,j);
        u_old(1,j) = u_old(2,j);
        
    end
    
    %V-cells
    for j=2:N
        
        v(1,j) = v(2,j);
        v_old(1,j) = v_old(2,j);
        
    end

%% North j=1 %%
    
    %T-cells
    for i=2:M-1
        
        eta(i,1) = eta(i,2);
        eta_old(i,1) = eta_old(i,2);
        
        H(i,1) = H(i,2);    
        H_old(i,1) = H_old(i,2);    
        
    end
    
    %U-cells
    for i=2:M
        
        u(i,1) = u(i,2);
        u_old(i,1) = u_old(i,2);
        
    end

    %V-cells
    for i=2:M-1
        
        v(i,1) = v(i,2);
        v_old(i,1) = v_old(i,2);

    end
    
%% East i=M+1(U) || M(T,V) %%
    
    %T-cells
    for j=2:N-1
        
        eta(M,j) = eta(M-1,j);
        eta_old(M,j) = eta_old(M-1,j);
        
        H(M,j) = H(M-1,j);    
        H_old(M,j) = H_old(M-1,j);    
                
    end
    
    %U-cells
    for j=2:N-1
        
        u(M+1,j) = u(M,j);
        u_old(M+1,j) = u_old(M,j);
        
    end
    
    %V-cells
    for j=2:N
                
        v(M,j) = v(M-1,j);
        v_old(M,j) = v_old(M-1,j);

    end

%% South j=N(T,U) OR N+1(V) %%

    %T-cells
    for i=2:M-1

        eta(i,N) = eta(i,N-1);
        eta_old(i,N) = eta_old(i,N-1);
        
        H(i,N) = H(i,N-1);    
        H_old(i,N) = H_old(i,N-1);    
        
    end

    %U-cells
    for i=2:M
                
        u(i,N) = u(i,N-1);
        u_old(i,N) = u_old(i,N-1);
        
    end
    
    %V-cells
    for i=2:M-1
                        
        v(i,N+1) = v(i,N);
        v_old(i,N+1) = v_old(i,N);

    end

%--------------------------------------------------------------------------
%SHEL SHallow-water numerical modEL
%Copyright (C) 2006,2009,2010. Guillaume Riflet, Instituto Superior Técnico
%da Universidade Técnica de Lisboa.
%--------------------------------------------------------------------------
    