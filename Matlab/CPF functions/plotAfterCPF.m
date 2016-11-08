function plotAfterCPF(systemName,caseName)
savename = sprintf('ResultsCPF/%s_%s_resultsCPF.mat',systemName,caseName);
resCPF1 = load(savename);
savename = sprintf('ResultsCPF/%s_%s_CPF.mat',systemName,caseName);
resCPF2 = load(savename);

mpc = openCase(systemName);
define_constants;
sysSettings = getSystemSettings(systemName,caseName);
rikt_vec = studyStochMod(systemName,caseName).';
rikt = rikt_vec(1,:);

nbIter = size(resCPF1.genQ,2);
nbLine = size(resCPF1.flowsAvg,1);
lineNames = cell(nbLine,1);
nbLoad = size(resCPF2.storeP,1);
loadNames = cell(nbLoad,1);
idxNZloads = find(mpc.bus(:,3) ~= 0);
idxLoad = idxNZloads(rikt~=0);
nbBus = size(resCPF2.storeV,1);
busNames = cell(nbBus,1);
nbGen = size(resCPF1.genQ,1);
genNames = cell(nbGen,1);

% Some flags to control the plotting
flagX = 'lambda';

% We want to single out the controllable generators and the slack.
slackbus = find(mpc.bus(:,2) == 3);
idxGenToPlot = [find(mpc.gen(:,1) == slackbus);sysSettings.indControl];

% Identify largest variations in V
% First, identify if we use the full CPF or stopped at nose point
idx_lower_part = find(resCPF2.storeK == -1,1,'first');
if isempty(idx_lower_part)
    dVlast = abs(resCPF2.storeV(:,end)-resCPF2.storeV(:,end-1));
else
    dVlast = abs(resCPF2.storeV(:,idx_lower_part)-resCPF2.storeV(:,idx_lower_part-1));
end
[~,sortedIdx_V] = sort(dVlast,'descend');
maxIdx_dV = sortedIdx_V(1:min(length(sortedIdx_V),5));

% Identify largest variations of Q
deltaQ = abs(resCPF1.genQ(:,end)-resCPF1.genQ(:,1));
[~,sortedIdx_Q] = sort(deltaQ,'descend');
maxIdx_dQ = sortedIdx_Q(1:min(length(sortedIdx_Q),5));
idxGenToPlot = sort(unique([idxGenToPlot;maxIdx_dQ]));

% Assemble predicted and corrected values
lambda_corr_pred = zeros(2*nbIter,1);
Vm_corr_pred = zeros(nbBus,2*nbIter);
lambda_corr_pred(1:2:end) = resCPF2.storeLambdaPred(1:2:end);
lambda_corr_pred(2:2:end) = resCPF2.storeLambda;
Vm_corr_pred(:,1:2:end) = resCPF2.storeVP(:,1:2:end);
Vm_corr_pred(:,2:2:end) = resCPF2.storeV;

% When did generators reach their limits?
iter_gen_hit_lim = find(resCPF2.storeLastHit);

% Construct matrix of values to plot for different parts of the PV curve
% (defined with respect to hitting the Qlim)
nb_hits = length(iter_gen_hit_lim);
colors_gen_hit_qlim = distinguishable_colors(nb_hits+1);

% How many different continuation parameters did we get?
map_cont_param_color = unique(resCPF2.storeK);
nb_cont_param = length(map_cont_param_color);
colors_cont_param = distinguishable_colors(nb_cont_param,[colors_gen_hit_qlim;1 1 1]);


for i = 1:nbLine
    linei = sprintf('Line %d -> %d.',mpc.branch(i,1),mpc.branch(i,2));
    lineNames{i} = linei;
end

for i = 1:nbLoad
    loadi = sprintf('Load at bus %d',i);
    loadNames{i} = loadi;
end

for i = 1:nbBus
    busi = sprintf('Bus %d',i);
    busNames{i} = busi;
end

for i = 1:nbGen
    nodegeni = mpc.gen(i,1);
    geni = sprintf('Gen at node %d',nodegeni);
    if mpc.bus(nodegeni,2) == 3
        geni = ['SLACK - ' geni];
    end
    genNames{i} = geni;
end

% fprintf('Maximum load values: %.2f\n',resCPF2.storeP(idxLoad,end));

% figure
% plot(resCPF2.storeLambda,imag(resCPF1.flowsAvg)*mpc.baseMVA);
% legend(lineNames);
% xlabel('Lambda');
% ylabel('Reactive power flows [MVA]');
% title([caseName ' - Reactive power flows [MVA]']);
% 
% figure
% plot(resCPF2.storeLambda,real(resCPF1.flowsAvg)*mpc.baseMVA);
% legend(lineNames);
% xlabel('Lambda');
% ylabel('Active power flows [MVA]');
% title([caseName ' - Active power flows [MVA]']);
% 

if strcmp(flagX,'lambda')
    xseries = resCPF2.storeLambda;
    xserlabel = 'lambda';
elseif strcmp(flagX,'nbiter')
    xseries = 1:nbIter;
    xserlabel = 'Nb of iterations';
end

% figure
% plot(xseries,resCPF2.storeStep);
% xlabel(xserlabel);
% ylabel('Step length');
% title([caseName ' - Step length']);
% 
% figure
% plot(xseries,resCPF2.storeK);
% xlabel(xserlabel);
% ylabel('Continuation parameter k');
% title([caseName ' - Continuation parameter']);

% figure
% plot(xseries,resCPF2.storeP(idxLoad,:)*mpc.baseMVA);
% xlabel(xserlabel);
% ylabel('Loads [MVA]');
% title([caseName ' - Loads [MVA]']);
% legend(loadNames(idxLoad));

figure
plot(xseries,resCPF2.storeLastHit,'o')
title('Generators reaching limits');

figure
hold on
plot(xseries,resCPF1.genQ*mpc.baseMVA);
genLeg = plot(xseries,resCPF1.genQ(idxGenToPlot,:)*mpc.baseMVA,'LineWidth',3);
plot([0 max(xseries)],repmat(mpc.gen(sysSettings.indControl,QMAX),1,2)'*mpc.baseMVA,'--');
xlabel(xserlabel);
ylabel('Reactive power generation [MVA]');
title([caseName ' - Qgen [MVA]']);
legend(genLeg,genNames(idxGenToPlot));

figure
hold on
plot(resCPF2.storeLambda,resCPF2.storeV);
plot(resCPF2.storeLambda,resCPF2.storeV(maxIdx_dV,:),'LineWidth',3);
xlabel('Lambda');
ylabel('Bus voltages [p.u.]');
title([caseName ' - Voltages']);
legend(busNames);

figure 
hold on
plot(lambda_corr_pred,Vm_corr_pred(sortedIdx_V(1),:),'r');
h_corr_all = plot(lambda_corr_pred(2:2:end),Vm_corr_pred(sortedIdx_V(1),2:2:end),'k','LineWidth',2);
% Plotting the different colors corresponding to different PV curves, as
% defined by generators having reached their limit.

if ~isempty(iter_gen_hit_lim)
    legend_array_qlim = zeros(nb_hits+1,1);
    legend_text_qlim = cell(nb_hits+1,1);
    for i = 0:nb_hits
        if i == 0
            iter_nbs = 1:iter_gen_hit_lim(1);
        elseif i == nb_hits
            iter_nbs = iter_gen_hit_lim(i):length(resCPF2.storeLambda);
        else
            iter_nbs = iter_gen_hit_lim(i):iter_gen_hit_lim(i+1);
        end
        hi = plot(resCPF2.storeLambda(iter_nbs),resCPF2.storeV(maxIdx_dV(1),iter_nbs),'Color',colors_gen_hit_qlim(i+1,:),'LineWidth',2);
        legend_array_qlim(i+1) = hi;
        if i == 0
            legend_text_qlim{i+1} = 'All below limits';
        else
            legend_text_qlim{i+1} = sprintf('Gen. %d at limit',resCPF2.storeLastHit(iter_gen_hit_lim(i)));
        end
    end
else
    legend_array_qlim = h_corr_all;
    legend_text_qlim{1} = 'PV curve';
end
% 
legend_array_cont_param = -10*ones(nb_cont_param,1);
legend_text_cont_param = cell(nb_cont_param,1);
for i = 1:2:length(lambda_corr_pred)
    % Plot predicted value, color depends on the continuation parameter for
    % that step.
    idx_color_pred = resCPF2.storeK((i+1)/2) == map_cont_param_color;
    color_pred = colors_cont_param(idx_color_pred,:);
    hi = plot(lambda_corr_pred(i),Vm_corr_pred(sortedIdx_V(1),i),...
        'o','MarkerSize',3,'MarkerEdgeColor',color_pred,'MarkerFaceColor',color_pred);
    if legend_array_cont_param(idx_color_pred) == -10
        % First time we encounter this parameter for plotting
        legend_array_cont_param(idx_color_pred) = hi;
        legend_text_cont_param{idx_color_pred} = sprintf('k = %d',...
            map_cont_param_color(idx_color_pred));
    end
end
xlabel('Loading \lambda');
ylabel(['Voltage magnitude at bus ' int2str(sortedIdx_V(1))]);
legend([legend_array_qlim;legend_array_cont_param],...
    {legend_text_qlim{:},legend_text_cont_param{:}});
title_txt = sprintf('Predictor-corrector process at bus %s, for dir = %s',...
    int2str(sortedIdx_V(1)),mat2str(resCPF2.dirCPF));
title(title_txt);
end