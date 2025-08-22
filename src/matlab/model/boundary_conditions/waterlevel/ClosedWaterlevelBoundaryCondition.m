classdef ClosedWaterlevelBoundaryCondition < IWaterlevelBoundaryCondition
    methods
        function state = apply(obj, state)
            % Implements closed boundary for waterlevel (eta)
            % Set boundary waterlevel to zero or land value
            state.eta([1 end],:) = 0;
            state.eta(:,[1 end]) = 0;
        end
    end
end
