classdef RadiationWaterlevelBoundaryCondition < IWaterlevelBoundaryCondition
    methods
        function state = apply(obj, state)
            % Implements radiation (Sommerfeld) boundary for waterlevel (eta)
            M = state.M; N = state.N;
            g = state.g;
            dx = state.dx;
            dt = state.dt;
            % East boundary
            state.eta(M,:) = state.eta(M,:) - dt/dx * sqrt(g * state.H(M,:)) .* (state.eta(M,:) - state.eta(M-1,:));
            % West boundary
            state.eta(1,:) = state.eta(1,:) - dt/dx * sqrt(g * state.H(1,:)) .* (state.eta(1,:) - state.eta(2,:));
            % Add North/South boundaries as needed
        end
    end
end
