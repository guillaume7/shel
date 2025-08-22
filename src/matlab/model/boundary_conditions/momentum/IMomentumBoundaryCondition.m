classdef (Abstract) IMomentumBoundaryCondition
    methods (Abstract)
        state = apply(obj, state);
    end
end
