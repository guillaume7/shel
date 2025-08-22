classdef (Abstract) ITracerBoundaryCondition
    methods (Abstract)
        state = apply(obj, state);
    end
end
