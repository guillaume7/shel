classdef NeumannMomentumBoundaryCondition < IMomentumBoundaryCondition
    methods
        function state = apply(obj, state)
            % Implements Neumann boundary for momentum
            % Example: copy interior values to boundary
            M = state.M; N = state.N;
            state.u(1,:) = state.u(2,:);
            state.u(M+1,:) = state.u(M,:);
            state.v(:,1) = state.v(:,2);
            state.v(:,N+1) = state.v(:,N);
        end
    end
end
