classdef NeumannTracerBoundaryCondition < ITracerBoundaryCondition
    methods
        function state = apply(obj, state)
            % Implements Neumann boundary for tracer
            M = state.M; N = state.N;
            state.Tr(1,:) = state.Tr(2,:);
            state.Tr(M,:) = state.Tr(M-1,:);
            state.Tr(:,1) = state.Tr(:,2);
            state.Tr(:,N) = state.Tr(:,N-1);
        end
    end
end
