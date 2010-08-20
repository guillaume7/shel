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