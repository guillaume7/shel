% Add solver subfolders to Matlab path
addpath('model/solvers');
addpath('model/solvers/time');
addpath('model/solvers/momentum');
addpath('model/solvers/tracer');
addpath('model/solvers/waterlevel');
addpath('model/boundary_conditions');
addpath('model/boundary_conditions/momentum');
addpath('model/boundary_conditions/tracer');
addpath('model/boundary_conditions/waterlevel');

% Modular entry point for SHEL
params = loadParams(); % You should implement loadParams to read config or GUI settings
state = ModelState(params);
state = InitialConditions.setup(state, params);
state = Grid.setup(state); % If needed, or merge with InitialConditions
solverType = params.solverType; % e.g., 'leapfrog'
numSteps = params.numSteps;

momentumSolver = SolverFactory.createMomentumSolver('upwind');
tracerSolver = SolverFactory.createTracerSolver('upwind');
solver = SolverFactory.createSolver('leapfrog', momentumSolver);

% Create boundary conditions using the factory
momentumBC = BoundaryConditionFactory.createMomentumBoundaryCondition(params.momentumBCType);
waterlevelBC = BoundaryConditionFactory.createWaterlevelBoundaryCondition(params.waterlevelBCType);
tracerBC = BoundaryConditionFactory.createTracerBoundaryCondition(params.tracerBCType);

for t = 1:numSteps
    state = solver.step(state);
    state = momentumBC.apply(state);
    state = waterlevelBC.apply(state);
    state = tracerSolver.step(state);
    state = tracerBC.apply(state);
    state = Diagnostics.compute(state, t);
    % Output, visualization, etc.
    OutputManager.saveOutput(state, t);
    OutputManager.plotModel(state);
    OutputManager.printStats(state, t);
end
