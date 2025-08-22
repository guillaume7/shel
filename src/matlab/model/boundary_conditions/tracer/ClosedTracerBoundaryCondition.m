classdef ClosedTracerBoundaryCondition < ITracerBoundaryCondition
    methods
        function state = apply(obj, state)
            % Implements closed boundary for tracer
            state.Tr([1 end],:) = 0;
            state.Tr(:,[1 end]) = 0;
        end
    end
end
