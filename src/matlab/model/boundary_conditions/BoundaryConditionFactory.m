classdef BoundaryConditionFactory
    methods (Static)
        function bc = createMomentumBoundaryCondition(type)
            switch type
                case 'closed'
                    bc = boundary_conditions.momentum.ClosedMomentumBoundaryCondition();
                case 'flather'
                    bc = boundary_conditions.momentum.FlatherMomentumBoundaryCondition();
                case 'neumann'
                    bc = boundary_conditions.momentum.NeumannMomentumBoundaryCondition();
                otherwise
                    error('Unknown momentum boundary condition type');
            end
        end
        function bc = createWaterlevelBoundaryCondition(type)
            switch type
                case 'closed'
                    bc = boundary_conditions.waterlevel.ClosedWaterlevelBoundaryCondition();
                case 'radiation'
                    bc = boundary_conditions.waterlevel.RadiationWaterlevelBoundaryCondition();
                case 'neumann'
                    bc = boundary_conditions.waterlevel.NeumannWaterlevelBoundaryCondition();
                otherwise
                    error('Unknown waterlevel boundary condition type');
            end
        end
        function bc = createTracerBoundaryCondition(type)
            switch type
                case 'closed'
                    bc = boundary_conditions.tracer.ClosedTracerBoundaryCondition();
                case 'radiation'
                    bc = boundary_conditions.tracer.RadiationTracerBoundaryCondition();
                case 'neumann'
                    bc = boundary_conditions.tracer.NeumannTracerBoundaryCondition();
                otherwise
                    error('Unknown tracer boundary condition type');
            end
        end
    end
end
