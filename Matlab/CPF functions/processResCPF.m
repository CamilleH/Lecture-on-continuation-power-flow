function processResCPF(systemName,caseName)
% This function can be run after ch_runCPF to process the results from CPF and
% compute additional values

savename = sprintf('ResultsCPF/%s_%s_CPF.mat',systemName,caseName);
resCPF = load(savename);
nbIterCPF = size(resCPF.storeP,2);
caseSettings = getSystemSettings(systemName,caseName);
mpc = openCase(systemName);
ngen = size(mpc.gen,1);
nline = size(mpc.branch,1);
nbus = size(mpc.bus,1);
Q_P = zeros(nbus,1);
Q_P(mpc.bus(:,3) ~= 0) = mpc.bus(mpc.bus(:,3) ~= 0,4)./mpc.bus(mpc.bus(:,3) ~= 0,3);
System = ch_system_init(mpc,caseSettings);
Ybus = System.Ybus_stat;
Yf = System.Yf;
Yt = System.Yt;

genP_int = mpc.gen(:,2);
if ~isempty(mpc.wind)
    Pwind = mpc.wind(:,2);
    injWind = zeros(nbus,1);
    injWind(mpc.wind(:,1)) = Pwind;
else
    injWind = 0;
end

busslack = mpc.bus(mpc.bus(:,2) == 3,1);

genQ = zeros(ngen,nbIterCPF);
injP = zeros(nbus,nbIterCPF);
injQ = zeros(nbus,nbIterCPF);
genPslack = zeros(nbIterCPF,1);
genQslack = zeros(nbIterCPF,1);
flowsFr = zeros(nline,nbIterCPF);
flowsTo = zeros(nline,nbIterCPF);

% For each iteration we want to compute additional values:
% - Production at slack
% - Reactive power productions at all generators
% - Flows on the lines

for i = 1:nbIterCPF
    Vbus = resCPF.storeVcomp(:,i);
    Pload = resCPF.storeP(:,i);
    Qload = Q_P.*Pload;
    Sload = Pload+1i*Qload;
    Str = Vbus.*conj(Ybus*Vbus);
    Snet = Str+Sload-injWind;
    
    % Net injection at each node
    injP(:,i) = real(Snet);
    injQ(:,i) = imag(Snet);
    
    % Flows on the line
    flowsFr(:,i) = Vbus(mpc.branch(:,1)).*conj(Yf*Vbus);
    flowsTo(:,i) = Vbus(mpc.branch(:,2)).*conj(Yt*Vbus);
    
    % Production of the slack bus
    genPslack(i) = real(Snet(busslack));
    
    % Q generation
    genQ(:,i) = injQ(mpc.gen(:,1),i);
end
flowsAvg = 1/2*(flowsFr-flowsTo);

savename = sprintf('ResultsCPF/%s_%s_resultsCPF.mat',systemName,caseName);
save(savename,'genQ','genPslack','flowsTo','flowsFr','flowsAvg');