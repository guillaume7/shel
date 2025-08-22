classdef NeumannWaterlevelBoundaryCondition < IWaterlevelBoundaryCondition
    methods
        function state = apply(obj, state)
            % Implements Neumann boundary for waterlevel (eta)
            M = state.M; N = state.N;
            state.eta(1,:) = state.eta(2,:);
            state.eta(M,:) = state.eta(M-1,:);
            state.eta(:,1) = state.eta(:,2);
            state.eta(:,N) = state.eta(:,N-1);
        end
    end
end
