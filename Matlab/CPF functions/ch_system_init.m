function System = ch_system_init(mpc,caseSettings)
%CH_SYSTEM_INIT Initializes the system values
%   Two possibilities: 
%     (i) Dynamical data for exciters and such is included
%     (ii) Dynamical data is not included
Settings_ch = ch_settings_init();
bifSettings = Settings_ch.bifSettings;
define_constants;

gen = mpc.gen;
bus = mpc.bus;
if isfield(mpc,'wind')
    wind = mpc.wind;
else
    wind = [];
end
branch = mpc.branch;

% Dynamical data if included

if isfield(mpc,'gen_dyn')
    gen_dyn = mpc.gen_dyn;
    ex = mpc.ex;
    System.dynData = 1;
else
    System.dynData = 0;
end
% System.dynData = 0;

nbus = size(bus,1);
indices.nbus = nbus;
indices.ng = size(gen,1);
indices.nwind = size(wind,1);
nlines = size(branch,1);
indices.nlines = nlines;

System.indices = indices;
System.gen = gen;
System.wind = wind;
System.baseMVA = mpc.baseMVA;

if System.dynData
    param.fs = mpc.fs;
    System.param = param;
    System.gen_dyn = gen_dyn;
    System.ex = ex;
end

%% Q_P ratio
indNZ = find(bus(:,PD));
Q_P = zeros(nbus,1);
Q_P(indNZ) = bus(indNZ,QD)./bus(indNZ,PD);
System.Q_P = Q_P;
System.incP = caseSettings.indLoads;
System.indWP = caseSettings.indWP;
System.indContSOPF = caseSettings.indControl; 
indUInLamb = [];
if bifSettings.inclContGen
    indUInLamb = [indUInLamb caseSettings.indControl];
end
System.indUInLamb = indUInLamb;

%% Y matrices

[Ybus_stat, Yf, Yt] = makeYbus(mpc);
if System.dynData
    Ybus_dyn = ch_calculateYbusDyn(gen,gen_dyn,Ybus_stat,System);
    System.Ybus_dyn = Ybus_dyn;
end

System.Ybus_stat = Ybus_stat;
System.Yf = Yf;
System.Yt = Yt;

%% Branch limit 
System.lineLim = branch(:,RATE_A)/mpc.baseMVA;

%% Connection matrix

% "From" matrix: nl*nb: one row for each branch, Clb = 1 if bus b is the
% fbus end of branch b
% "To" matrix: nl*nb: one row for each branch, Clb = 1 if bus b is the tbus
% end of branch b

indB = 1:nlines;
indF = branch(:,F_BUS);
indT = branch(:,T_BUS);
valuesC = ones(nlines,1);

Cf = sparse(indB,indF,valuesC,nlines,nbus);
Ct = sparse(indB,indT,valuesC,nlines,nbus);

System.Cf = Cf;
System.Ct = Ct;
System.indBranchF = indF;
System.indBranchT = indT;
System.branch = branch;
