classdef (Abstract) IWaterlevelBoundaryCondition
    methods (Abstract)
        state = apply(obj, state);
    end
end
