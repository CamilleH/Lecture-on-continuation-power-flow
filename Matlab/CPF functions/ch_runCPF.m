function resultsCPF = ch_runCPF(systemName,caseName,c,dirCPF)

% system: system name
% u: power generation of the generators
% c = contingency number

%% In this file, we perform many CPF with different values of the generations u

% System settings
%eval([system '_settings;']);

% Load increase direction
if nargin == 2
    c=0;
    caseSettings = getSystemSettings(systemName,caseName);
    dir = caseSettings.rikt_vec_all(:,1);
elseif nargin == 3
    %     dir_vec = studyStochMod(systemName,caseName).';
    %     dir = dir_vec(1,:);
    caseSettings = getSystemSettings(systemName,caseName);
    dir = caseSettings.rikt_vec_all(:,1);
elseif nargin == 4
    dir = dirCPF;
else
    error('Wrong number of inputs');
end

% Running the CPFs
optionsCPF.dirCPF = dir/norm(dir);
% When comparing the CPF with optimization problem, make sure to set
% chooseStartPoint to 0
optionsCPF.chooseStartPoint = 0; 
% optionsCPF.chooseStartPoint = 1;
optionsCPF.startPoint = 0;

try
    resultsCPF = ch_CPF_Dyn(systemName,caseName,c,optionsCPF);
catch
    keyboard;
end
% try
%     resultsCorr = ch_processAfterCPF(systemName,caseName,c,optionsCPF,resultsCPF);
% catch exception
%     keyboard;
%     fprintf('\n The process after the CPF failed.\n');
%     fprintf('The error was:\n%s\n\n',exception.getReport());
% end

% Extracting results
[nbIter,storeP,storeV,storeVcomp,storeVP,storeEig,storeLambda,storeLambdaPred,...
    storeLastHit,storeK,storeStep] = ...
    v2struct(resultsCPF,{'fieldNames','nbIter','storeP','storeV','storeVcomp','storeVP','storeEig','storeLambda','storeLambdaPred',...
    'storeLastHit','storeK','storeStep'});

% Cleaning up the results from the CPF and saving
storeP = storeP(:,1:nbIter-1);
storeV = storeV(:,1:nbIter-1);
storeVcomp = storeVcomp(:,1:nbIter-1);
storeVP = storeVP(:,1:2*(nbIter-1));
storeEig = storeEig(:,1:(nbIter-1));
storeLambda = storeLambda(1:nbIter-1);
storeLambdaPred = storeLambdaPred(1:2*(nbIter-1));
storeLastHit = storeLastHit(1:nbIter-1);
storeK = storeK(1:nbIter-1);
storeStep = storeStep(1:nbIter-1);

toSave = 1;
if toSave
    savename = sprintf('ResultsCPF/%s_%s_CPF.mat',systemName,caseName);
    save(savename,'dirCPF','storeP','storeV','storeVcomp','storeVP','storeEig','storeLambda','storeLambdaPred','storeLastHit','storeK','storeStep');
end
end