classdef (Abstract) IMomentumSolver
    methods (Abstract)
        coreRHSu = computeCoreRHSu(obj, state);
        coreRHSv = computeCoreRHSv(obj, state);
    end
    methods
        function RHSu = computeRHSu(obj, state)
            coreRHSu = obj.computeCoreRHSu(state);
            drag_u = obj.bottomDragU(state);
            wind_u = obj.windStressU(state);
            RHSu = coreRHSu + wind_u - drag_u;
        end
        function RHSv = computeRHSv(obj, state)
            coreRHSv = obj.computeCoreRHSv(state);
            drag_v = obj.bottomDragV(state);
            wind_v = obj.windStressV(state);
            RHSv = coreRHSv + wind_v - drag_v;
        end
        function drag_u = bottomDragU(obj, state)
            if isfield(state, 'Cd') && ~isempty(state.Cd)
                H_u = MomentumUtilities.twoaverage_u(state.H, state.M, state.N);
                drag_u = MomentumForcing.bottom_drag(state.u, H_u, state.Cd);
            else
                drag_u = zeros(size(state.u));
            end
        end
        function drag_v = bottomDragV(obj, state)
            if isfield(state, 'Cd') && ~isempty(state.Cd)
                H_v = MomentumUtilities.twoaverage_u(state.H', state.N, state.M)';
                drag_v = MomentumForcing.bottom_drag(state.v, H_v, state.Cd);
            else
                drag_v = zeros(size(state.v));
            end
        end
        function wind_u = windStressU(obj, state)
            if isfield(state, 'uwind') && isfield(state, 'rho0') && isfield(state, 'rho_air') && ~isempty(state.uwind)
                wind_u = MomentumForcing.wind_stress(state.u, state.rho0, state.rho_air, state.uwind, 'u');
            else
                wind_u = zeros(size(state.u));
            end
        end
        function wind_v = windStressV(obj, state)
            if isfield(state, 'vwind') && isfield(state, 'rho0') && isfield(state, 'rho_air') && ~isempty(state.vwind)
                wind_v = MomentumForcing.wind_stress(state.v, state.rho0, state.rho_air, state.vwind, 'v');
            else
                wind_v = zeros(size(state.v));
            end
        end
    end
end
