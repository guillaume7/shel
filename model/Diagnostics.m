classdef Diagnostics
    methods (Static)
        function state = compute(state, l)
            % Compute diagnostics (energy, vorticity, etc.)
            M = state.M; N = state.N; dx = state.dx; dy = state.dy;
            dA = dx * dy;
            % Kinetic energy
            u_t = .5 * (state.u(1:M,:) + state.u(2:M+1,:));
            v_t = .5 * (state.v(:,1:N) + state.v(:,2:N+1));
            Hnoland = state.H;
            Hnoland(Hnoland < -1) = 0;
            state.ke = .5 * state.rho0 * dA * (u_t.^2 + v_t.^2) .* Hnoland .* state.mask;
            state.iKe(l) = sum(sum(state.ke));
            % Potential energy
            state.pe = .5 * state.rho0 * state.g * dA .* state.eta.^2 .* state.mask;
            state.iPe(l) = sum(sum(state.pe));
            % Eddy kinetic energy (simplified)
            state.gradux = (state.u(2:M+1,:) - state.u(1:M,:)) / dx;
            state.gradvy = (state.v(:,2:N+1) - state.v(:,1:N)) / dy;
            state.graduy = zeros(M,N); state.gradvx = zeros(M,N);
            state.graduy(2:M,2:N-2) = (state.u(2:M,3:N-1) - state.u(2:M,2:N-2)) / dy;
            state.gradvx(2:M-2,2:N) = (state.v(3:M-1,2:N) - state.v(2:M-2,2:N)) / dx;
            state.eke = state.eke + state.rho0 * dA * state.K * 2 * state.dt * (state.gradux.^2 + state.gradvy.^2 + state.graduy.^2 + state.gradvx.^2) ...
                .* Hnoland .* state.mask;
            state.iEke(l) = sum(sum(state.eke));
            % Volume
            state.volume(l) = dA * sum(sum(state.eta .* state.mask));
            % Momentum
            state.iMomentumU(l) = dA * sum(sum(u_t .* Hnoland,2));
            state.iMomentumV(l) = dA * sum(sum(v_t .* Hnoland,1));
            % Add more diagnostics as needed
        end
    end
end
