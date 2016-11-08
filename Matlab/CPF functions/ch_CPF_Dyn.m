function resultsCPF = ch_CPF_Dyn(systemName,caseName,c,optionsCPF)
%% File for the dynamic CPF
define_constants;

if c==0
    sysToLoad = systemName;
else
    sysToLoad = sprintf('%s_cont%d',systemName,c);
end

% load the file
mpc = openCase(sysToLoad);
% mpc.gen(:,QMAX) = 9999;
caseSettings = getSystemSettings(systemName,caseName);

% Matpower constants
define_constants;

% define some global settings
Settings_ch = ch_settings_init();
System = ch_system_init(mpc,caseSettings);

% Get the settings
cpfSettings = Settings_ch.cpfSettings;
maxIter = cpfSettings.maxIter;
stepLambda = cpfSettings.stepLambda;
stepV = cpfSettings.stepV;
verbose = cpfSettings.verbose;

% Get the indices
indices = System.indices;
nbus = indices.nbus;
ngen = indices.ng;
lineLim = System.lineLim;

% Indices of the Efp in y
indX_Efp = (3*ngen+1):(4*ngen);

% Normalize the load and generation values from the Matpower case
bus = mpc.bus;
gen = mpc.gen;
if isfield(mpc,'wind')
    wind = mpc.wind;
else
    wind = [];
end
baseMVA = mpc.baseMVA;

% Direction
indNZ = find(bus(:,PD));
Q_P = zeros(nbus,1);
Q_P(indNZ) = bus(indNZ,QD)./bus(indNZ,PD);
dirP = zeros(nbus,1);
dirP(caseSettings.indLoads) = optionsCPF.dirCPF;
dirQ = dirP.*Q_P;

%% generator info
% on = find(gen(:, GEN_STATUS) > 0);      %% which generators are on?
% gbus = gen(on, GEN_BUS);                %% what buses are they at?

% Initialiazing the values
V = bus(:,VM).*exp(1i*bus(:,VA));
% [ref, pv, pq] = bustypes(bus, gen);
% vcb = ones(size(V));           %% create mask of voltage-controlled buses
% vcb(pq) = 0;                    %% exclude PQ buses
% idx_gen_pv = find(vcb(gbus));            %% in-service gens at v-c buses
% V(gbus(idx_gen_pv)) = gen(on(idx_gen_pv), VG) ./ abs(V(gbus(idx_gen_pv))).* V(gbus(idx_gen_pv));
% V = ones(size(bus,1),1);
if System.dynData
    delta = zeros(size(gen,1),1);
    omega = zeros(size(gen,1),1);
    Eqp = ones(size(gen,1),1);
    Ef = ones(size(gen,1),1);
    x = [delta;omega;Eqp;Ef];
else
    x = [];
end

% Choosing another starting point if wanted
if optionsCPF.chooseStartPoint > 0
    bus(caseSettings.indLoads,PD) = optionsCPF.startPoint;
    bus(:,QD) = bus(:,PD).*Q_P;
end

if optionsCPF.chooseStartPoint > 1 || optionsCPF.chooseStartPoint == -1
    % Then, we also set the new generation at some of the gens
    if ~isempty(caseSettings.indControl) && ~isempty(optionsCPF.startPointU)
        gen(caseSettings.indControl,PG) = optionsCPF.startPointU;
    end
    if ~isempty(wind) && ~isempty(optionsCPF.startPointWP)
        wind(:,PG) = optionsCPF.startPointWP;
    end
end
P0 = bus(:,PD);
Q0 = bus(:,QD);
P_cur = bus(:,PD);
Q_cur = bus(:,QD);
if ~isempty(wind)
    Pwind0 = wind(:,2);
else
    Pwind0 = [];
end

lambda = 0;
nbIter = 1;
type = 0;
%success = 1;

% Preparing the arrays for storing the results
storeP = zeros(nbus,maxIter);
storeV = zeros(nbus,maxIter);
storeVcomp = zeros(nbus,maxIter);
storeVP = zeros(nbus,2*maxIter);
storeEig = zeros(ngen*4-2,maxIter);
storeLambda = zeros(maxIter,1);
storeLambdaPred = zeros(2*maxIter,1);
storeLastHit = zeros(maxIter,1);
storeSkipIter = zeros(maxIter,1);
storeStep = zeros(maxIter,1);
storeK = zeros(maxIter,1);

% Initializing with a PF
allSettings{1} = Settings_ch;
allSettings{2} = System;
if System.dynData
    resultsPFini = ch_PF_Dyn(x,V,gen,bus,wind,1,allSettings);
    [successIni,x,V,bus,gen_a,gen_b,last_hit,f_x,f_y,g_x,g_y] = v2struct(resultsPFini,{'fieldNames','success','x','V','bus','gen_a','gen_b','last_hit','f_x','f_y','g_x','g_y'});
else
    resultsPFini = ch_PF_Dyn(x,V,gen,bus,wind,1,allSettings);
    [successIni,x,V,bus,gen_a,gen_b,last_hit] = v2struct(resultsPFini,{'fieldNames','success','x','V','bus','gen_a','gen_b','last_hit'});
    % Update the gen matrix from the results
    gen(gen_b,8) = 0; % deactivate all PQ generators
end
redStepPred = 0;
flagForcek = 0;
stepPred = stepLambda;
indTransLim = [];
if successIni
    k = ch_choosePara(System,Settings_ch,x,V,dirP,dirQ,gen_a,gen_b,bus);
    success = 1;
    while nbIter <= maxIter && success
        if verbose > 0
            fprintf(1,'--------- CPF Iteration %d ------------\n',nbIter);
        end
        
        % Sets pv and pq
        pv = bus(:,2) == 2;
        pq = bus(:,2) == 1;
        npv = sum(pv);
        npq = sum(pq);
        
        % Predictor
        if redStepPred
            stepPred = stepPred/10;
        end

        [xpred,Vpred,lambdapred,t] = ch_predictDyn(System,x,V,lambda,gen_a,gen_b,dirP,dirQ,stepPred,bus,k);
        storeVP(:,2*(nbIter-1)+1) = abs(Vpred);
        storeLambdaPred(2*(nbIter-1)+1) = lambdapred;
        
        % Updating the power
        P_cur = bus(:,PD);
        Q_cur = bus(:,QD);
        bus(:,PD) = P0+lambdapred*dirP;
        bus(:,QD) = Q0+lambdapred*dirQ;

        % Storing the tangent vector for later use
        if System.dynData
            Efp_b = xpred-x;
            Efp_b = Efp_b(indX_Efp);
            tangentCont = [xpred-x;Efp_b;angle(Vpred-V);abs(Vpred-V)];
        else
            tangentCont = [angle(Vpred-V);abs(Vpred-V)];
            tangentCont = t(1:end-1);
        end
        
        % Choosing the parameter k
        % in the tangent vector, we hav first terms for Va and then Vm
        if System.dynData
            indVmInt = 4*ngen-2+nbus+(1:nbus+1);
        else
            indVmInt = npv+npq+1:length(t);
        end
        
        dV_dLambda = t(indVmInt);
        indVm = find(pq);
        [tsorted,idxSorted] = sort(abs(dV_dLambda),'descend');
        if ~flagForcek
            % if flagForcek == 1, we have forced using some specific k so
            % we don't have to change it here.
            k = idxSorted(1);
            if k == length(dV_dLambda)
                k = 0;
                dV_max = tsorted(2)/tsorted(1); % Largest voltage decrease
                if dV_max > cpfSettings.slopeLim1
                    % Case where we force k to be one of the Vms because
                    % the slope is large enough
                    k = idxSorted(2);
                end
                if nbIter > 2 && lambda < storeLambda(nbIter-2)
                    % case where we follow the lower part of the nose curve
                    % i.e. decreasing value of lambda
                    k = -1;
                end
            end
            if k > 0 && ~System.dynData
                k = indVm(k);
            end
        end
        storeK(nbIter) = k;
        
        % corrector
        % storing the gen a and b
        gen_a_cur = gen_a;
        gen_b_cur = gen_b;
        
        last_hit = [];

        try
            if System.dynData
                resultsPFcorr = ch_PF_Dyn(xpred,Vpred,gen,bus,wind,0,allSettings,dirP,dirQ,k,gen_a,gen_b);
                [success,x_corr,V_corr,bus_corr,deltaLambda_loc,gen_a,gen_b,last_hit,JacDyn,f_x,f_y,g_x,g_y] = v2struct(resultsPFcorr,{'fieldNames',...
                    'success','x','V','bus','deltaLambda_loc','gen_a','gen_b','last_hit','JacDyn','f_x','f_y','g_x','g_y'});
            else
                resultsPFcorr = ch_PF_Dyn(xpred,Vpred,gen,bus,wind,0,allSettings,dirP,dirQ,k,gen_a,gen_b);
                [success,x_corr,V_corr,bus_corr,deltaLambda_loc,gen_a,gen_b,last_hit,JacDyn] = v2struct(resultsPFcorr,{'fieldNames',...
                    'success','x','V','bus','deltaLambda_loc','gen_a','gen_b','last_hit','JacDyn'});
            end
        catch
            keyboard;
        end
        
        % Set the parameter for next step
        if ~success && k == 0  %|| k==-1)
            % case where we don't find a solution and the parameter is still
            % lambda, we try to force the parameter to be the critical voltage,
            % even if the derivative of the latter is rather large. Then we
            % retry.
            % We use the values from the previous step
            gen_a = gen_a_cur;
            gen_b = gen_b_cur;
            bus(:,PD) = P_cur;
            bus(:,QD) = Q_cur;
            if System.dynData
                k = idxSorted(2);
            else
                k = indVm(idxSorted(2));
            end
%             k = ch_choosePara(x,V,dirP,dirQ,gen_a,gen_b,bus,1);
            flagForcek = 1;
            success = 1;
            redStepPred = 1;
            continue;
        elseif ~success && stepPred > cpfSettings.stepMin
            gen_a = gen_a_cur;
            gen_b = gen_b_cur;
            bus(:,PD) = P_cur;
            bus(:,QD) = Q_cur;
            success = 1;
            redStepPred = 1;
            continue;
        elseif ~success && stepPred < cpfSettings.stepMin
            break;
%         else
%             [k,slopeVmin] = ch_choosePara(x_corr,V_corr,dirP,dirQ,gen_a,gen_b,bus);
        end

        if length(last_hit) == 1
            storeLastHit(nbIter) = last_hit;
        end
              
        if length(last_hit) > 1
            % We reduce the time steps if two generators have reached their
            % limits at the same time
            redStepPred = 1;
            bus(:,PD) = P_cur;
            bus(:,QD) = Q_cur;
            gen_a = gen_a_cur;
            gen_b = gen_b_cur;
            storeSkipIter(nbIter) = 1;
            nbIter = nbIter + 1;
            success = 1;
            continue;
        elseif ~System.dynData
            % Adjusting the gen based on results from power flows
            gen(last_hit,8) = 0;
            bus(:,2) = bus_corr(:,2); % updating pv, pq types
        end
        
        flagForcek = 0;
        
        if k>0
            lambda_corr = lambdapred+deltaLambda_loc;
        else
            lambda_corr = lambdapred;
        end
        
        % Check if the CPF is finished
        % Is lambda negative?
        if lambda_corr < 0 || k==-1 && lambda_corr < 0.5
            fprintf(1,'Lambda is negative. END OF CPF.\n');
            break;
        end       
        
        % Computing the eigenvalues
        D = eig(full(JacDyn));
        
        % Extracting the one with maximum real part
        if ~cpfSettings.hopf
            Dcheck = D(imag(D)==0);
        else
            Dcheck = D;
        end
        
        if System.dynData
            [Dmax_real,ind_max] = max(real(Dcheck));
            if isempty(Dmax_real)
                % case when we don't consider Hopf bifurcations
                % (cpfSettings.hopf = 0) but all eigenvalues are purely
                % imaginary.
                Dmax_real = -99;
            end
            Dmax_real = Dmax_real(1);
            condBif = Dmax_real > 0;
        else
            % We look at the eigenvalues of JacPF
            [Dmin,ind_min] = min(abs(Dcheck));
            Dmin = Dmin(1);
            condBif = Dmin < 0.01;
        end
        
        % Attempt: compute the tangent at the corrected point to see
        % whether it is now pointing toward smaller lambdas, in which case
        % it is an SLL
        if ~isempty(last_hit)
            [xpred2,Vpred2,lambdapred2] = ch_predictDyn(System,x_corr,V_corr,lambda_corr,gen_a,gen_b,dirP,dirQ,stepPred,bus,k);
            dirLambda = sign(lambdapred2-lambda_corr);
            condBif = condBif | dirLambda == -1;
        end

        if ~cpfSettings.allTheWay && (condBif || ~success || lambda_corr < lambda)
            break;
        end
        
        % Checking whether some lines have reached their limits
        [Sf,St] = ch_calculateLineFlows(V_corr,System);
        Pf = real(Sf);
        Pt = real(St);
        Pmean = abs(Pf-Pt)/2;
        indTransLim = find(abs(Pf)>lineLim,1);
        if verbose > 1
            PmeanMW = Pmean*baseMVA;
            lineLimMW = lineLim*baseMVA;
            fprintf(1,'-- flows on the lines \n--');
            %         fprintf(1,'Active flow on branch %d: %.2f MW and limit %2.f MW:\n',...
            %             [1:length(Pmean);PmeanMW.';lineLimMW.']);
            fprintf(1,'Active flow on branch %d: %.2f MW and limit %2.f MW:\n',...
                [1:length(Pf);baseMVA*Pf.';lineLimMW.']);
        end
        if ~cpfSettings.allTheWay && cpfSettings.flowLimits && ~isempty(indTransLim)
            type = 4;
            break;
        end
        
        % Step adaptation
        if cpfSettings.adaptStep
            cpf_error = norm([[angle(V_corr(pq));abs(V_corr(pv | pq))]-[angle(Vpred(pq));abs(Vpred(pv |pq))];deltaLambda_loc],inf);
            resultsPFcorr.nbItersPerStep(resultsPFcorr.nbItersPerStep == 0) = [];
            if cpf_error < cpfSettings.errorTol
                %% Increase stepsize
                stepPred = stepPred*cpfSettings.errorTol/cpf_error;
                if stepPred > cpfSettings.stepMax
                    stepPred = cpfSettings.stepMax;
                end
            elseif isempty(last_hit) && resultsPFcorr.nbItersPerStep(end) > 5
                %% decrese stepsize
                stepPred = stepPred*cpfSettings.errorTol/cpf_error;
                if stepPred < cpfSettings.stepMin
                    stepPred = cpfSettings.stepMin;
                end
            end
        else
            if k == 0 || k == -1
                stepPred = stepLambda;
            else
                stepPred = stepV;
            end
        end
        redStepPred = 0;
        
        % Updating
        x = x_corr;
        V = V_corr;
        lambda = lambda_corr;
        
        bus_bak = bus;
        bus(:,PD) = P0+lambda_corr*dirP;
        bus(:,QD) = Q0+lambda_corr*dirQ;
        
        if ~System.dynData
            deltaQ = computeGenQ(System,V,bus);
            indLim_all = findExceedGenQ(deltaQ,gen,1:size(gen,1));
            % Compare with gen_b
            genbLim = ismember(find(gen_b),indLim_all);
            if sum(genbLim) ~= length(genbLim)
                % We come here if a generator is now able to produce enough
                % reactive power
                keyboard;
            end
        end
        storeLambda(nbIter) = lambda_corr;
        storeV(:,nbIter) = abs(V_corr);
        storeVcomp(:,nbIter) = V_corr;
        storeVP(:,2*(nbIter-1)+2) = abs(V_corr);
        storeLambdaPred(2*(nbIter-1)+2) = lambda_corr;
        storeP(:,nbIter) = bus(:,PD);
        storeStep(nbIter) = stepPred;
        
        storeEig(1:length(D),nbIter) = D;
        nbIter = nbIter+1;
    end
    
    %% Going back one step to get the last stable point before bifurcation
    bus(:,PD) = P_cur;
    bus(:,QD) = Q_cur;
end

if verbose > 0
    fprintf(1,'CPF completed in %d / %d.\n',nbIter,maxIter);
end 

if nbIter == maxIter
    keyboard;
end

% Store results
if ~System.dynData
    f_x=[];f_y=[];g_x=[];g_y=[];
end
resultsCPF = v2struct(nbIter,storeP,storeV,storeVcomp,storeVP,storeEig,storeLambda,storeLambdaPred,...
    storeLastHit,storeK,storeStep,JacDyn,type,gen_a,gen_b,xpred,Vpred,gen,wind,bus,x,V,tangentCont,...
    last_hit,indTransLim,f_x,f_y,g_x,g_y,lambda,lambda_corr,successIni);