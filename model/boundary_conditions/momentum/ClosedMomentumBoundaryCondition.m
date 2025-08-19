classdef ClosedMomentumBoundaryCondition < IMomentumBoundaryCondition
    methods
        function state = apply(obj, state)
            % Implements closed boundary for momentum (u, v)
            % Set boundary velocities to zero
            state.u([1 end],:) = 0;
            state.u(:,[1 end]) = 0;
            state.v([1 end],:) = 0;
            state.v(:,[1 end]) = 0;
        end
    end
end
