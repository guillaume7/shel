classdef ContinuitySolver
    methods (Static)
        function RHSeta = computeRHSeta(state)
            % Compute the right-hand side of the continuity equation (waterlevel)
            % Implements divergence of fluxes on the Arakawa C-grid
            M = state.M; N = state.N;
            dx = state.dx; dy = state.dy;
            mask = state.mask;
            u = state.u; v = state.v;
            RHSeta = zeros(M,N);
            % Compute flux divergence (centered differences)
            RHSeta(2:M-1,2:N-1) = - mask(2:M-1,2:N-1) .* ( ...
                (u(3:M,2:N-1) - u(2:M-1,2:N-1)) / dx ...
                + (v(2:M-1,3:N) - v(2:M-1,2:N-1)) / dy ...
            );
            % Add atmospheric pressure deviation forcing
            if isfield(state, 'p_surf') && ~isempty(state.p_surf)
                dH = PressureForcing.pressure_deviation(state.p_surf, state.rho0, state.g);
                RHSeta(2:M-1,2:N-1) = RHSeta(2:M-1,2:N-1) + dH(2:M-1,2:N-1);
            end
            % Open boundary treatment can be added here (e.g., radiation, Flather)
            % For now, boundaries are left as zeros (can be expanded)
        end
    end
end
