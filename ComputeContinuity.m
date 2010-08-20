function RHSeta = ComputeContinuity(H, eta, u, v, mask)
%function ComputeContinuity
%Shallow-water equations solver

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

% global H; 
% global eta; 
% global u; 
% global v; 
% global mask; 
% global RHSeta;
global dx; 
global dy; 

[M,N] = size(eta);
RHSeta = zeros(M,N); 

%Compute along U then along V
RHSeta(2:M-1,2:N-1) = - mask(2:M-1,2:N-1) .* ( ...
    ComputeContinuityU(H,u,mask,dx) + ComputeContinuityU(H',v',mask',dy)' ...
);

function RHS = ComputeContinuityU(H,u,mask,dx)
[M,N] = size(H);
RHS = ( ... 
            mask(3:M,2:N-1) .* .5 ...
            .* ( H(3:M,2:N-1) + H(2:M-1,2:N-1) ) .* u(3:M,2:N-1) ... %i+1, i, i+1 
            - mask(1:M-2,2:N-1) .* .5 ...
            .* ( H(2:M-1,2:N-1) + H(1:M-2,2:N-1) ) .* u(2:M-1,2:N-1) ... %i, i-1, i 
        )/dx;