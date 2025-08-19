% Add solver subfolders to Matlab path
addpath('model/initial_conditions');
addpath('model/initial_conditions/bottom');
addpath('model/initial_conditions/momentum');
addpath('model/initial_conditions/tracer');
addpath('model/initial_conditions/waterlevel');
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
M = state.M; N = state.N;

% Initialize bottom bathymetry
state.d = BottomInitialFactory.create(params.bottomICType, M, N, params);
% Initialize waterlevel
state.eta = WaterlevelInitialFactory.create(params.waterlevelICType, M, N, params);
state.eta_old = state.eta;
state.eta_new = state.eta;
state.H = state.eta + state.d;
state.H_old = state.H;
state.H_new = state.H;
% Initialize momentum
state.u = MomentumInitialFactory.create('constant', M+1, N, state.u0);
state.u_old = state.u;
state.u_new = state.u;
state.u_a = state.u;
state.v = MomentumInitialFactory.create('constant', M, N+1, state.v0);
state.v_old = state.v;
state.v_new = state.v;
state.v_a = state.v;
% Initialize tracer
state.Tr = TracerInitialFactory.create(params.tracerICType, M, N, params);

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
    % Compute diagnostics
    state = LocalDiagnostics.computeAll(state);
    state = GlobalDiagnostics.computeAll(state, t);
    % Output, visualization, etc.
    OutputManager.saveOutput(state, t);
    OutputManager.plotModel(state);
    OutputManager.printStats(state, t);
end
