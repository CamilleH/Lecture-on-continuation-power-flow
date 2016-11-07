% Script to compare CPF and repetitive load flows
close all;clc;clear;
define_constants;
mpc0 = loadcase('case39');
nonzero_loads = mpc0.bus(:,PD)~=0;
tanphi = mpc0.bus(nonzero_loads,QD)./mpc0.bus(nonzero_loads,PD);

%% CPF
mpopt = mpoption('verbose', 2);
mpopt = mpoption(mpopt, 'cpf.stop_at', 'NOSE');
mpopt = mpoption(mpopt,'cpf.plot.level',2,'cpf.parameterization',2,'cpf.adapt_step',1);
mpct = mpc0;
mpct.bus(nonzero_loads,PD) = mpct.bus(nonzero_loads,PD)+500; 
mpct.bus(nonzero_loads,QD) = tanphi.*mpct.bus(nonzero_loads,PD); 
results_cpf = runcpf(mpc0,mpct,mpopt);

%% Repeated load flows
deltaP = 1;
possible = 1;
mpc = mpc0;
mpopt = mpoption('out.all',0);
pvfig = figure;
while possible
    % We increase the loads by deltaP MW at constant power factors
    mpc.bus(nonzero_loads,PD) = mpc.bus(nonzero_loads,PD)+deltaP;
    mpc.bus(nonzero_loads,QD) =tanphi.*mpc.bus(nonzero_loads,PD);
    results = runpf(mpc,mpopt);
    possible = results.success;
    % Plot
    if possible
        figure(pvfig);
        hold on
        plot(results.bus(12,PD),results.bus(12,VM),'*');
    end
end

%% Comparison
[mpc.bus(:,PD) results_cpf.bus(:,PD)]
