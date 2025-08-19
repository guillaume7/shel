classdef UpwindMomentumSolver < IMomentumSolver
    methods
        function coreRHSu = computeCoreRHSu(obj, state)
            M = state.M; N = state.N;
            dx = state.dx; g = state.g; f = state.f;
            mask_u = state.mask_u;
            u = state.u; v = state.v;
            H = state.H; H_old = state.H_old;
            v_old = state.v_old; v_u = zeros(size(u));
            H_u = zeros(size(u)); H_old_u = zeros(size(u));
            % Averaging for staggered grid
            H_u(2:M,2:N-1) = MomentumUtilities.twoaverage_u(H,M,N);
            H_old_u(2:M,2:N-1) = MomentumUtilities.twoaverage_u(H_old,M,N);
            v_old_u(2:M,2:N-1) = MomentumUtilities.fouraverage_u(v_old(:,2:N),M,N-2);
            v_u(2:M,2:N-1) = MomentumUtilities.fouraverage_u(v(:,2:N),M,N-2);
            coreRHSu = zeros(size(u));
            coreRHSu(2:M,2:N-1) = mask_u(2:M,2:N-1) .* ( ...
                - g * (H_u(2:M,2:N-1) - H_u(1:M-1,2:N-1)) / dx ...
                + f * .5 * (v_u(2:M,2:N-1) + v_old_u(2:M,2:N-1)) ...
            );
        end
        function coreRHSv = computeCoreRHSv(obj, state)
            M = state.M; N = state.N;
            dy = state.dy; g = state.g; f = state.f;
            mask_v = state.mask_v;
            u = state.u; v = state.v;
            H = state.H; H_old = state.H_old;
            u_old = state.u_old; u_v = zeros(size(v));
            H_v = zeros(size(v)); H_old_v = zeros(size(v));
            H_v(2:M-1,2:N) = MomentumUtilities.twoaverage_u(H',N,M)';
            H_old_v(2:M-1,2:N) = MomentumUtilities.twoaverage_u(H_old',N,M)';
            u_old_v(2:M-1,2:N) = MomentumUtilities.fouraverage_u(u_old(2:M,:),N-2,M)';
            u_v(2:M-1,2:N) = MomentumUtilities.fouraverage_u(u(2:M,:),N-2,M)';
            coreRHSv = zeros(size(v));
            coreRHSv(2:M-1,2:N) = mask_v(2:M-1,2:N) .* ( ...
                - g * (H_v(2:M-1,2:N) - H_v(2:M-1,1:N-1)) / dy ...
                - f * .5 * (u_v(2:M-1,2:N) + u_old_v(2:M-1,2:N)) ...
            );
        end
    end
end
