function [RHSu, RHSv] = ComputeUV2D_old
%% function newu = ComputeV2D
%Shallow-water equations solver
%
% Centred differences
%
%NOTE: Any viscous and frictional term
%must be evaluated in tindex(i-1,M)e l-1.
%
% switch v <-> u; j <-> i and dy <-> dx relative to ComputeU2D
% switch v_old <-> u_old
% take special care for the coriolis force term
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

global H_old; %to use in dissipation terms 
global H; 
global eta; 
global u; 
global u_old; %to use in dissipation terms 
global v; 
global v_old; 
global dx; 
global dy; 
global M; 
global N; 
global mask_u; 
global mask_v; 
%global mask; 
global g; 
global f; 
global K; 
global uwind; 
global vwind; 
global wind; 
global bottom; 
global pressure; 
global coriolis; 
global wall;

global rho0;
global rho_air;
global lb; 
global karman; 

wall = false; 

%% U spatial scheme

RHSu = zeros(size(u));
H_u = zeros(size(u));
H_old_u = zeros(size(u));
v_old_u = zeros(size(u));
v_u = zeros(size(u));
H_u(2:M,2:N-1) = twoaverage_u(H,M,N);
H_old_u(2:M,2:N-1) = twoaverage_u(H_old,M,N);
v_old_u(2:M,2:N-1) = fouraverage_u(v_old(:,2:N),M, N-2); 
v_u(2:M,2:N-1) = fouraverage_u(v(:,2:N),M, N-2); 

RHSu(2:M,2:N-1) = mask_u(2:M,2:N-1) .* ( ... 
    ...%advection along U ... 
        - 0.2500 / dx * ( ... %i-1, i 
        mask_u(3:M+1,2:N-1) .* ( u(3:M+1,2:N-1) + u(2:M,2:N-1) ).^2 .* H(2:M,2:N-1) ... %i+1, i, i+1, i, i 
      - mask_u(1:M-1,2:N-1) .* ( u(2:M,2:N-1) + u(1:M-1,2:N-1) ).^2 .* H(1:M-1,2:N-1) ... %i, i-1, i, i-1, i-1 
     )  ... 
    ...%advection along V 
     - 0.0625 / dy * ( ... %i-1, i 
        ( ... 
            mask_u(2:M,3:N) .* ( H(2:M,3:N) + H(2:M,2:N-1) + H(1:M-1,3:N) + H(1:M-1,2:N-1) ) ... % j+1 ; j ; i-1,j+1 ; i-1 
          .* ( v(2:M,3:N) + v(1:M-1,3:N) ) .* ( u(2:M,3:N) + u(2:M,2:N-1) ) ...% j+1 ; i-1,j+1 ; j+1 ; i 
          - mask_u(2:M,1:N-2) .* ( H(2:M,2:N-1) + H(2:M,1:N-2) + H(1:M-1,2:N-1) + H(1:M-1,1:N-2) ) ... % i ; j-1 ; i-1 ; i-1,j-1 
          .* ( v(2:M,2:N-1) + v(1:M-1,2:N-1) ) .* ( u(2:M,2:N-1) + u(2:M,1:N-2) ) ...% i ; i-1,j-1 ; i ; j-1 
        ) ... 
     ) ... 
    ...%Coriolis Force 
    + coriolis * ( ... 
        f * v_u(2:M,2:N-1) .* H_u(2:M,2:N-1) ... % i ; i-1 . 
    )... 
    ...%Pressure gradient Force 
    + pressure * ( ... 
       - g * H_u(2:M,2:N-1) .* ( eta(2:M,2:N-1) - eta(1:M-1,2:N-1) ) / dx ...  
                               ... % i ; i-1 ;                 i ; i-1 
    )... 
    ...%Viscosity along U 
      + ( ... % i ; i-1 
          + K * mask_u(3:M+1,2:N-1) .* H_old(2:M,2:N-1) .* ( u_old(3:M+1,2:N-1) - u_old(2:M,2:N-1) ) / dx ... % i ; i+1 ; i ; i 
          - K * mask_u(1:M-1,2:N-1) .* H_old(1:M-1,2:N-1) .* ( u_old(2:M,2:N-1) - u_old(1:M-1,2:N-1) ) / dx ... % i-1 ; i ; i-1; i-1 
      ) / dx ... % i ; i-1 
    ...%Viscosity along V 
      + ( ... % i ; i-1 ; i ; i-1 . 
          + K / dy * mask_u(2:M,3:N) .* H_old(2:M,2:N-1) .* ( u_old(2:M,3:N) - u_old(2:M,2:N-1) ) ... % (j;j+1;i-1,j+1;i-1);j+1 ; j 
          - K / dy * mask_u(2:M,1:N-2) .* H_old(2:M,1:N-2) .* ( u_old(2:M,2:N-1) - u_old(2:M,1:N-2) ) ... % (j-1;j;i-1;i-1,j-1);j ; j-1 
      ) / dy ...
      );
 
    %bottom stress. Look out not to stress out more than there is...
    if bottom || wall 
        stress = zeros(size(RHSu(2:M,2:N-1))); 
        if bottom 
            stress = bottomstress( u_old(2:M,2:N-1), v_old_u(2:M,2:N-1), ...
                        H_old_u(2:M,2:N-1), ...
                        lb, karman ); 
        end 
        if wall 
            stress = stress + bottomstress(u_old(2:M,2:N-1), 0, dy, lb, karman) .* ( ... 
                (1 - mask(1:M-1,3:N)) .* H_old(1:M-1,3:N) + ... 
                (1 - mask(2:M,3:N)) .* H_old(2:M,3:N) + ... 
                (1 - mask(1:M-1,1:N-2)) .* H_old(1:M-1,1:N-2) + ... 
                (1 - mask(2:M,1:N-2)) .* H_old(2:M,1:N-2) ... 
            ) * 0.5 / dy; 
        end 
        ind = find( (RHSu(2:M,2:N-1) - mask_u(2:M,2:N-1) .* stress) ./ RHSu(2:M,2:N-1) > 0 );
        RHSu(ind) = RHSu(ind) - mask_u(ind) .* stress(ind); 
        ind = find( (RHSu(2:M,2:N-1) - mask_u(2:M,2:N-1) .* stress) ./ RHSu(2:M,2:N-1) <= 0 );
        RHSu(ind) = 0.; 
    end
    
%Wind stress
if wind 
    RHSu(2:M,2:N-1) = RHSu(2:M,2:N-1) ...
            + mask_u(2:M,2:N-1) .* windstress(...
                uwind * ones(size(u(2:M,2:N-1))), ...
                vwind * ones(size(u(2:M,2:N-1))), ...
                rho0, rho_air ...
            );
end
    
%% V spatial scheme

%Local variables for hydrodynamics 
RHSv = zeros(size(v));
H_v = zeros(size(v));
H_old_v = zeros(size(v));
u_old_v = zeros(size(v));
u_v = zeros(size(v));
H_v(2:M-1,2:N) = twoaverage_v(H,M,N); 
H_old_v(2:M-1,2:N) = twoaverage_v(H,M,N); 
u_old_v(2:M-1,2:N) = fouraverage_v(u_old(2:M,:), M-2, N); 
u_v(2:M-1,2:N) = fouraverage_v(u(2:M,:), M-2, N);

RHSv(2:M-1,2:N) = mask_v(2:M-1,2:N) .* ( ... 
    ...%advection along V ... 
        - 0.2500 / dy * ( ... %j-1, j 
        mask_v(2:M-1,3:N+1) .* ( v(2:M-1,3:N+1) + v(2:M-1,2:N) ).^2 .* H(2:M-1,2:N) ... %j+1, j, j+1, j, j 
      - mask_v(2:M-1,1:N-1) .* ( v(2:M-1,2:N) + v(2:M-1,1:N-1) ).^2 .* H(2:M-1,1:N-1) ... %j, j-1, j, j-1, j-1 
     )  ... 
    ...%advection along U 
     - 0.0625 / dx * ( ... %j-1, j 
        ( ... 
            mask_v(3:M,2:N) .* ( H(3:M,2:N) + H(2:M-1,2:N) + H(3:M,1:N-1) + H(2:M-1,1:N-1) ) ... % i+1 ; i ; i+1, j-1 ; j-1 
          .* ( u(3:M,2:N) + u(3:M,1:N-1) ) .* ( v(3:M,2:N) + v(2:M-1,2:N) ) ...% i+1 ; i+1, j-1 ; i+1 ; j 
          -  mask_v(1:M-2,2:N) .*( H(2:M-1,2:N) + H(1:M-2,2:N) + H(2:M-1,1:N-1) + H(1:M-2,1:N-1) ) ... % j ; i-1 ; j-1 ; i-1, j-1 
          .* ( u(2:M-1,2:N) + u(2:M-1,1:N-1) ) .* ( v(2:M-1,2:N) + v(1:M-2,2:N) ) ...% j ; i-1, j-1 ; j ; i-1 
        ) ... 
     ) ... 
    ...%Coriolis Force 
    + coriolis * ( ... 
    - f * u_v(2:M-1,2:N) .* H_v(2:M-1,2:N) ... 
    )... 
    ...%Pressure gradient Force 
    + pressure * ( ... 
       - g * H_v(2:M-1,2:N) .* ( eta(2:M-1,2:N) - eta(2:M-1,1:N-1) ) / dy ... 
                               ... % j ; j-1 ;                 j ; j-1 
    )... 
    ...%Viscosity along U 
      + ( ... % j ; j-1 
          + K * mask_v(2:M-1,3:N+1) .* H_old(2:M-1,2:N) ...
          .* ( v_old(2:M-1,3:N+1) - v_old(2:M-1,2:N) ) / dy ... % j ; j+1 ; j ; j 
          - K * mask_v(2:M-1,1:N-1) .* H_old(2:M-1,1:N-1) ...
          .* ( v_old(2:M-1,2:N) - v_old(2:M-1,1:N-1) ) / dy ... % j-1 ; j ; j-1; j-1 
      ) / dy ... % j ; j-1 
    ...%Viscosity along V 
      + ( ... % j ; j-1 ; j ; j-1 . 
          + K / dx * mask_v(3:M,2:N) .* H_old(2:M-1,2:N) ...
          .* ( v_old(3:M,2:N) - v_old(2:M-1,2:N) ) ... % (i;i+1;i+1, j-1;j-1);i+1 ; i 
          - K / dx * mask_v(1:M-2,2:N) .* H_old(1:M-2,2:N) ...
          .* ( v_old(2:M-1,2:N) - v_old(1:M-2,2:N) ) ... % (i-1;i;j-1;i-1, j-1);i ; i-1 
      ) / dx ...
     );
 
    %bottom stress. Look out not to stress out more than there is... 
    if bottom || wall 
        stress = zeros(size(RHSv(2:M-1,2:N))); 
        if bottom 
            stress = bottomstress( v_old(2:M-1,2:N), u_old_v(2:M-1,2:N), ...
                                H_old_v(2:M-1,2:N), lb, karman ); 
        end 
        if wall 
            stress = stress + bottomstress( v_old(2:M-1,j), 0, dx, ...
                            lb, karman ) .* ( ... 
                (1 - mask(i+1,j-1)) .* H_old(i+1,j-1) + ... 
                (1 - mask(i+1,  j)) .* H_old(i+1,  j) + ... 
                (1 - mask(i-1,j-1)) .* H_old(i-1,j-1) + ... 
                (1 - mask(i-1,  j)) .* H_old(i-1,  j) ... 
            ) * 0.5 / dx; 
        end 
        ind = find( (RHSv(2:M-1,2:N) -  mask_v(2:M-1,2:N) .* stress) ./ RHSv(2:M-1,2:N) > 0 );
        RHSv(ind) = RHSv(ind) - mask_v(2:M-1,2:N) .* stress(ind); 
        ind = find( (RHSv(2:M-1,2:N) -  mask_v(2:M-1,2:N) .* stress) ./ RHSv(2:M-1,2:N) <= 0 );
        RHSv(ind) = 0.; 
    end 
    
%Wind stress
if wind 
    RHSv(2:M-1,2:N) = RHSv(2:M-1,2:N) ...
            + mask_v(2:M-1,2:N) .* windstress(...
                vwind * ones(size(v(2:M-1,2:N))), ...
                uwind * ones(size(v(2:M-1,2:N))), ...
                rho0, rho_air ...
            );
end

function av = fouraverage_u(v, M, N) 
%% function av = fouraverage_u(v, mask, M, N)

    av =  .25 * ( ... 
            v(2:M,1:N) + ... 
            v(2:M,2:N+1) + ... 
            v(1:M-1,1:N) + ... 
            v(1:M-1,2:N+1) ... 
           ); 
  
function av = fouraverage_v(u, M, N) 
%% function av = fouraverage_v(u, mask) 
    av = .25 * ( ... 
            u(1:M,2:N) + ... 
            u(2:M+1,2:N) + ... 
            u(1:M,1:N-1) + ... 
            u(2:M+1,1:N-1) ... 
           ); 
       
function av = twoaverage_u( H, M, N)
    av = .5 * ( ...
            H(2:M,2:N-1) + H(1:M-1,2:N-1) ...
         );
        
function av = twoaverage_v( H, M, N)
    av = .5 * ( ...
            H(2:M-1,2:N) + H(2:M-1,1:N-1) ...
         );

function tb_L = bottomstress( u_L, v_L, H_L, lb, karman) 
%% Bottom stress        
%function tb_L = bottomstress( u_L, v_L, H_L ) 
% 
%H_L is the total depth height 
     
    cd_L = karman / log( ( 0.5 .* H_L + lb ) / lb ) ; 
    % Checked-out Blumberg and Mellor 1987 
    % Double-checked with Pietrzak 2002 
    cd_L = max(cd_L .* cd_L, 0.0025); 
    tb_L = cd_L .* sqrt ( u_L .* u_L + v_L .* v_L ) .* u_L; 
     
function u_L = windstress(u_L, v_L, rho0, rho_air) 
%% Wind stress            
%function tv_L = windstress(v_L) 

    vmod_L = sqrt( u_L .* u_L + v_L .* v_L );

    if vmod_L < 6.
        Cd_L = -.26097 .* vmod_L +  2.4316;
    else
        Cd_L = .064986 .* vmod_L + .44053;
    end

    Cd_L = Cd_L * .001;

    u_L = Cd_L .* rho_air / rho0 .* vmod_L .* u_L;
