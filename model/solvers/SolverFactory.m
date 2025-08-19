classdef SolverFactory
    methods (Static)
        function momentumSolver = createMomentumSolver(type)
            switch type
                case 'upwind'
                    momentumSolver = momentum_solvers.UpwindMomentumSolver();
                % Add more cases for other momentum solvers
                otherwise
                    error('Unknown momentum solver type');
            end
        end
        function tracerSolver = createTracerSolver(type)
            switch type
                case 'upwind'
                    tracerSolver = tracer_solvers.UpwindTracerSolver();
                % Add more cases for other tracer solvers
                otherwise
                    error('Unknown tracer solver type');
            end
        end
        function solver = createSolver(type, momentumSolver)
            switch type
                case 'leapfrog'
                    solver = solvers.LeapfrogSolver(momentumSolver);
                % Add more cases for other solvers
                otherwise
                    error('Unknown solver type');
            end
        end
    end
end
