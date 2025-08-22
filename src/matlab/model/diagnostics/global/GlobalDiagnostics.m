classdef GlobalDiagnostics
    methods(Static)
        function state = computeAll(state, t)
            % Compute all global diagnostics and store in state
            dA = state.dx * state.dy;
            Hnoland = state.H;
            Hnoland(Hnoland < -1) = 0;
            state.volume(t) = dA * sum(state.eta(:));
            state.iMomentumU(t) = dA * sum(sum(0.5 * (state.u(1:end-1,:) + state.u(2:end,:)) .* Hnoland));
            state.iMomentumV(t) = dA * sum(sum(0.5 * (state.v(:,1:end-1) + state.v(:,2:end)) .* Hnoland));
            state.iKe(t) = sum(state.ke(:)) * dA;
            state.iPe(t) = sum(state.pe(:)) * dA;
            state.iEke(t) = sum(state.eke(:)) * dA;
            state.iVort(t) = sum(state.curl_t(:)) * dA;
            state.iEnst(t) = sum(state.enstrophy_t(:)) * dA;
            state.iSqShear(t) = sum(state.sqshearrate_t(:)) * dA;
            state.iSqStrech(t) = sum(state.sqstrechrate_t(:)) * dA;
            state.iSqStrain(t) = sum(state.sqstrain_t(:)) * dA;
            state.iOWeiss(t) = sum(state.okuboweiss_t(:)) * dA;
            state.vtime(t) = state.time;
        end
    end
end
