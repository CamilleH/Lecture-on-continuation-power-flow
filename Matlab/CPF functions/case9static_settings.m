function caseSettings = case9static_settings(caseName)
% This file contains the settings to be used in ch_main together with
% Karystianos' system

% open case
mpc = openCase('case9static');
define_constants;
% The indices of the controllable generators 
indControl = [2 3];
indGenCont = ismember(mpc.gen(:,1),indControl);
% indWP = 2;

% List of contingency number
contingencies = 0;

% The loads at which the consumption is increased
if strcmp(caseName,'load5')
    indLoads = 5;
elseif strcmp(caseName,'load6')
    indLoads = 6;
elseif strcmp(caseName,'load8')
    indLoads = 8;
elseif strcmp(caseName,'loads579')
    indLoads = [5 7 9];
end

%% Defining the direction of load increase.
nloads = length(indLoads);
rikt_vec_all = eye(nloads);
if nloads>1
    rikt_vec_all = [rikt_vec_all;ones(1,nloads)];
end
if strcmp(caseName,'loads579')
%     rikt_vec_all = [4.535;2.404;1.729]-mpc.bus([5;7;9],PD);
%     rikt_vec_all = [2.582;1.702;2.321]-mpc.bus([5;6;8],PD);
%     rikt_vec_all = rikt_vec_all/norm(rikt_vec_all);
    rikt_vec_all = [1;1;1];
end

%% Starting points

% Defining the starting point as the mean of the forecast
stochMod = case9static_StochModel(caseName,1);
startWP = stochMod.muBeta.';
startLoad = stochMod.muNorm;
startU = mpc.gen(indGenCont,2).';
% startU = 0;

if ~isempty(mpc.wind)
    indWP = mpc.wind(:,1);
else
    indWP = [];
end

caseSettings.indControl = find(indGenCont);
caseSettings.indWP = indWP;
caseSettings.indLoads = indLoads;
caseSettings.caseName = caseName;
caseSettings.startLoad = startLoad;
caseSettings.startU = startU;
caseSettings.startWP = startWP;
caseSettings.rikt_vec_all = rikt_vec_all;
caseSettings.contingencies = contingencies;