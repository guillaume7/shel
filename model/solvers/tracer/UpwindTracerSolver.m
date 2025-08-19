classdef UpwindTracerSolver < ITracerSolver
    methods
        function state = step(obj, state)
            M = state.M; N = state.N; dx = state.dx; dy = state.dy; dt = state.dt;
            K = state.K;
            RHSt_x = Utilities.ComputeSpaceT_UP_U(state.H, state.u_a, state.mask_u, state.Tr, K, dx);
            RHSt_y = Utilities.ComputeSpaceT_UP_U(state.H', state.v_a', state.mask_v', state.Tr', K, dy)';
            RHSt = RHSt_x + RHSt_y;
            minRHSt = min(min(RHSt));
            ind = find(RHSt == minRHSt);
            dttr = min( floor( abs(- state.H(ind) .* state.Tr(ind) ./ minRHSt) / dt ), 1000) * dt;
            state.Tr_new = state.mask .* (state.H .* state.Tr + dttr * RHSt) ./ state.H_new;
            state.Tr = state.Tr_new;
        end
    end
end
