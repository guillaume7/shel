classdef LeapfrogSolver < ISolver
    properties
        momentumSolver
    end
    methods
        function obj = LeapfrogSolver()
            obj.momentumSolver = MomentumSolver();
        end
        function state = step(obj, state)
            % Main leapfrog + Asselin-Roberts filter logic
            % Compute continuity (waterlevel update)
            RHSeta = ContinuitySolver.computeRHSeta(state);
            state.eta_new = state.eta_old + 2 * state.dt * RHSeta;
            % Compute new H
            state.H_new = state.eta_new + state.d;
            % Time stepping for u and v (vectorized, simplified)
            state.u_new = state.u_old + 2 * state.dt * ...
                obj.momentumSolver.computeRHSu(state);
            state.v_new = state.v_old + 2 * state.dt * ...
                obj.momentumSolver.computeRHSv(state);
            % Boundary conditions
            state = BoundaryConditions.apply(state);
            % Asselin-Roberts filter
            gama = state.gama;
            state.eta = state.eta + gama * (state.eta_old - 2 * state.eta + state.eta_new);
            state.u = state.u + gama * (state.u_old - 2 * state.u + state.u_new);
            state.v = state.v + gama * (state.v_old - 2 * state.v + state.v_new);
            % Update old/new
            state.eta_old = state.eta;
            state.eta = state.eta_new;
            state.H_old = state.H;
            state.H = state.H_new;
            state.u_old = state.u;
            state.u = state.u_new;
            state.v_old = state.v;
            state.v = state.v_new;
        end
    end
end
