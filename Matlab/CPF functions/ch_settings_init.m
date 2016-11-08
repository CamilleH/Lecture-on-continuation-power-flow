
function Settings_ch = ch_settings_init()

%% Settings for the power flow
pfSettings.tol = 1e-8;
pfSettings.maxIter = 20;
pfSettings.verbose = 0;
pfSettings.qlimFlag = 1; % Should we take into account the reactive power limits ?
pfSettings.qlimMode = 1; % How do we handle the reactive power limits? 1: the largest becomes a PQ bus

%% Settings for the continuation power flow
cpfSettings.maxIter = 10000;
cpfSettings.stepLambda = 0.1;
cpfSettings.normPredict = 0; % Whether we should normalize the step
cpfSettings.stepV = 0.01;
cpfSettings.allTheWay = 0;
cpfSettings.adaptStep = 0;
cpfSettings.stepMax = 0.2;
cpfSettings.stepMin = 5e-5;
cpfSettings.errorTol = 1e-3;
cpfSettings.verbose = 0;
cpfSettings.flagCorr = 1; % correction method (1: parametrization, 2: orthogonal)
cpfSettings.slopeLim1 = 0.2; % threshold for the V/lambda slope, switch 1-2
cpfSettings.slopeLim2 = 0.2; % threshold for the V/lambda slope, switch 2-3
cpfSettings.startZero = 0; % Whether we start from zero load power or from the initial load flow
cpfSettings.plotPV = 0;
cpfSettings.hopf = 0; % Do we consider Hopf bifurcations?
cpfSettings.printPV = 0; % Print the PV-curve
cpfSettings.flowLimits = 0; % consider the flow limits or not

%% Settings for calculating the distance to the closest bifurcation point
bifSettings.maxIter = 10000;
bifSettings.adaptStep = 0;
bifSettings.heavyBall = 1;
bifSettings.impFunc = 'quad'; %'dens' 'quad'
bifSettings.verbose = 1;
bifSettings.tol = 1e-4;
bifSettings.stepSurf = 0.01;
bifSettings.stepSurf2 = 0.005;
bifSettings.miniStep = 0.001; % minimum value of the step length
bifSettings.improvThre = 5e-3;
bifSettings.relImprovThre = 1e-3;
bifSettings.modeChangeStep = 2; % 1 or 2 depending on if we adapt the step length depending on the criteria or the distance
bifSettings.redStepSize = 0; % say whether we try to reapproach a surface with smaller step sizes
bifSettings.allTheWay = 1; % Continue on the next surface when reaching an edge or not
bifSettings.showSearch = 0; % indicates whether we plot the search along the path
% bifSettings.closestU = 1; % find the closest approximation point in (u,zeta)
bifSettings.inclContGen = 0; % Include the controllable generators in lambda
bifSettings.inclWP = 1; % Include the wind farms in lambda
% Mode: explore or comparison
% 1. explore: in all directions, find the path to the closest point, even if
% it is on a surface already explored.
% 2. comparison: keep just the surfaces that have not already been explored
bifSettings.modeSimu = 2;
bifSettings.includeSLL = 1;

%% Settings for the overall method
resSettings.verboseRes = 0;
resSettings.saveResults = 1;
resSettings.printPath = 0;
resSettings.plotNewSurf = 0;

%% Settings for MC simulations
mcSettings.storeXopt = 0;
mcSettings.storeSamples = 0;

%% Storing in the global variable
Settings_ch.pfSettings = pfSettings;
Settings_ch.cpfSettings = cpfSettings;
Settings_ch.bifSettings = bifSettings;
Settings_ch.resSettings = resSettings;
Settings_ch.mcSettings = mcSettings;
