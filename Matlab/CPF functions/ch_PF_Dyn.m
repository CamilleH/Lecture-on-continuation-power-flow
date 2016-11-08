function results = ch_PF_Dyn(x0,V0,gen,bus,wind,mode,allSettings,dirP,dirQ,k,gen_a,gen_b)
% Run power flows for the generators with one-axis model and AVR
% In x are: - delta for all generators but the slack
%           - omega for all gen (not anymore)
%           - Eqp for all gen
%           - Efp for the non limited exciters (set a)
%           - Efp for the limited exciters (set b)
%
% In V are: - Va for all normal buses
%           - Vm for all normal buses
% define_constants;

Settings_ch = allSettings{1};
System = allSettings{2};
pfSettings = Settings_ch.pfSettings;
tol = pfSettings.tol;
maxIter = pfSettings.maxIter;
verbose = pfSettings.verbose;

indices = System.indices;
if System.dynData
    ex = System.ex;
end

if nargin == 7 % case normal power flows (not corrector step of CPF)
    k = 0;
else
    P0 = bus(:,3);
    Q0 = bus(:,4);
end
if ~isempty(wind)
    Pwind = wind(:,2);
else
    Pwind = [];
end
deltaLambda_loc = 0;

%% Indices
% global indices
ng = indices.ng;
nbus = indices.nbus;
nwind = indices.nwind;

if System.dynData
    indZ_U = 4*ng-2+nbus+(1:nbus);
    indF_P = 4*ng-2+(1:nbus);
    indF_Q = 4*ng-2+nbus+(1:nbus);
    n_z = length(x0)+2*nbus-2;
else
    indZ_U = nbus+(1:nbus);
    n_z = 2*nbus;
    indF_P = 1:nbus;
    indF_Q = nbus+(1:nbus);
end

% Numerical indices of the sets a and b (AVRs that have not, and have,
% reached their limit, respectively);
if nargin == 7
    gen_a = ones(ng,1); % all exciters are assumed to be within the limits in the beginning
    gen_b = ~gen_a;
end

gen_a_num = find(gen_a);
gen_b_num = find(gen_b);
ngen_a = length(gen_a_num);

%% Parameters
if System.dynData
    Ef_lim = ex(:,9);
end

% Calculate JacDyn or not (it is computationally demanding so if we can
% avoid doing it, we avoid doing it)
flagNoJacDyn = mode == 2;

%% Beginning of the power flow
finished = 0;
success = 1;
last_hit = [];
%[~,JacDyn,~] = ch_calculateJacDyn(x0,V0,gen_a,gen_b,System);
nbStep = 0;
nbItersPerStep = zeros(10,1);

while ~finished
    nbIter = 1;
    Va = angle(V0);
    Vm = abs(V0);
    x = x0;
    V = V0;
    nbStep = nbStep+1;
    % Power flow
    
    % Put here to initialize them in case deltaF < threshold from the
    % beginning
%     if System.dynData
%         [JacPF,JacDyn,f_x,f_y,g_x,g_y] = ch_calculateJacDyn(x,V,gen_a,gen_b,System);
%     else
%         [JacPF,JacDyn] = ch_calculateJacDyn(x,V,gen_a,gen_b,System,bus);
%     end
    
    while nbIter <= maxIter
        
        if verbose > 0
            fprintf(1,'Iteration number %d.\n',nbIter);
        end

        % Below, if we have only static data (only power equations), F is
        % empty
        [F,G] = ch_calculateFGDyn(x,V,bus,gen_a,gen_b,gen(:,2),Pwind,System);

        if k > 0
            Vk = abs(V(k));
            deltaV = Vk-abs(V0(k));
            deltaF = [F;G;deltaV];
        else
            deltaF = [F;G];
        end
        if norm(deltaF,inf) < tol
            if ~exist('JacDyn','var') && nbIter == 1
                % Here we are already close enough from the solution at the
                % first iteration. 
                if System.dynData
                    [JacPF,JacDyn,f_x,f_y,g_x,g_y] = ch_calculateJacDyn(x,V,gen_a,gen_b,System,bus,flagNoJacDyn);
                else
                    [JacPF,JacDyn] = ch_calculateJacDyn(x,V,gen_a,gen_b,System,bus,flagNoJacDyn);
                end
            end
            break;
        end
        
        % Jacobian, warning: delta and omega for the slack have been left
        % out of the Jacobian. Thus, the increment will have one dimension
        % less than the number of variables
        % If only power flow equations are of interest, JacDyn is empty
        if System.dynData
            [JacPF,JacDyn,f_x,f_y,g_x,g_y] = ch_calculateJacDyn(x,V,gen_a,gen_b,System,bus,flagNoJacDyn);
        else
            [JacPF,JacDyn] = ch_calculateJacDyn(x,V,gen_a,gen_b,System,bus,flagNoJacDyn);
        end
        if k > 0
            if System.dynData
                ne = n_z+1;
                indenz = indZ_U(k);
            else
                pv = bus(:,2) == 2;
                pq = bus(:,2) == 1;
                npv = nnz(pv);
                npq = nnz(pq);
                ne =npv+2*npq+1;
                % The nonzero element in e is the k:th voltage magnitude
                % among all PQ buses
                indPqInVm = find(pq);
                indenz = npv+npq+find(k==indPqInVm);
            end
            e = zeros(1,ne);
            e(indenz) = 1;
            F_lambda = zeros(n_z,1);
            F_lambda(indF_P) = dirP;
            F_lambda(indF_Q) = dirQ;
            if ~System.dynData
                slack = pv==0&pq==0;
                indRemove = [slack;~pq];
                F_lambda(indRemove) = [];
            end
            A = [JacPF F_lambda;
                e];
        else
            A = JacPF;
        end
%         if rcond(full(A)) < 1e-15
%             keyboard;
%         end
        inc = -A\deltaF;
        if System.dynData
            deltaX = inc(1:(4*ng-2));
            deltaY = inc((4*ng-1):(4*ng+2*nbus-2));
        else
            nbInc = length(inc);
            deltaX = [];
            deltaY = inc;
        end
        if k > 0
            deltaLambda = inc(end);
            deltaLambda_loc = deltaLambda_loc+deltaLambda;
            if ~System.dynData
                % case where we consider dlambda as well and we have to
                % remove it from deltaY since we put all inc in deltaY
                % above
                deltaY(end) = [];
            end
        end
        
        if System.dynData
            % Warning: Ef occupy the last ng elements of x but are not sorted
            % according to the bus number. Instead, the Ef in gen_a comes first and
            % the Ef in gen_b comes second
            deltaEf_a = deltaX((3*ng-1):(3*ng+ngen_a-2));
            deltaEf_b = deltaX((3*ng+ngen_a-1):(4*ng-2));
            
            % Updating delta, omega and Eqp
            x(2:ng) = x(2:ng) + deltaX(1:(ng-1));
            x((ng+2):(2*ng)) = x((ng+2):(2*ng)) + deltaX(ng:(2*ng-2)); %
            %         OMEGA: not anymore
            x((2*ng+1):(3*ng)) = x((2*ng+1):(3*ng)) + deltaX((2*ng-1):(3*ng-2));
            % Updating Efp_a
            x(3*ng+gen_a_num) = x(3*ng+gen_a_num) + deltaEf_a;
            % Updating Efp_b
            x(3*ng+gen_b_num) = x(3*ng+gen_b_num) + deltaEf_b;
        end
        % Updating the voltages
        if System.dynData
            Va = Va + deltaY(1:nbus);
            Vm = Vm + deltaY((nbus+1):(2*nbus));
        else
            % We update Va in non PV or slack nodes and Vm in non slack
            % nodes
            % First entries in deltaY are for Va and then comes Vm
            pv = bus(:,2) == 2;
            pq = bus(:,2) == 1;
            npvpq = sum(pv)+sum(pq);
            Va(pv|pq) = Va(pv|pq) + deltaY(1:npvpq);
            Vm(pq) = Vm(pq) + deltaY((npvpq+1):end);
        end
        V = Vm.*exp(1i*Va);
        
        if k>0
            bus(:,3) = bus(:,3) + deltaLambda*dirP;
            bus(:,4) = bus(:,4) + deltaLambda*dirQ;
        end
        
        nbIter = nbIter+1;
    end
    
    iterArrayLength = length(nbItersPerStep);
    if iterArrayLength < nbStep
        nbItersPerStep(iterArrayLength+1:iterArrayLength+10) = 0;
    end
    nbItersPerStep(nbStep) = nbIter;
    
    if nbIter > maxIter
        % We don't even check the reactive power limits if the power flow
        % has not converged, because the divergence might give strange
        % values.
        success = 0;
        break;
    end
    
    % Computing the reactive power productions
    
    % Checking whether some exciters have reached their limits
    if System.dynData
        Ef = x((3*ng+1):(4*ng));
        diffEf = abs(Ef) > Ef_lim;
        indLim = find(diffEf);
    else
        % first compute the reactive power production of the generators
        Ybus_stat = System.Ybus_stat;
        Str = V.*conj(Ybus_stat*V);
        SL = bus(:,3) + 1i*bus(:,4);
        deltaQ = imag(SL) + imag(Str);
        gen_slack = bus(gen(:,1),2) == 3;
        % We don't check the reactive power at the slack bus or at the PV
        % buses that have become PQ buses.
        on = find(gen_a & ~gen_slack);
%         if isempty(on)
%             keyboard;
%         end
        Qg = deltaQ(gen(on,1));
        Qmax = gen(on,4);
        diffQ = Qg > Qmax;
        indLim = find(diffQ);
    end
    
    if isempty(indLim) || isempty(find(gen_a == true, 1))
        break;
    else
        if mode == 0 && length(indLim) > 1
            % mode == 1 corresponds to a normal PF (not inside CPF). In
            % that case, we consider even cases where several generators
            % reached their limits.
            last_hit = indLim;
            break;
        end
        if System.dynData
            [~,last_hit] = max(abs(Ef) - Ef_lim);
            gen_b(last_hit) = 1;
            gen_a = ~gen_b;
            gen_a_num = find(gen_a == true);
            gen_b_num = find(gen_b == true);
            ngen_a = length(gen_a_num);
        else
            [~,last_hit] = max(Qg - Qmax);
            % Changing the node from PV to PQ
            indLim2 = on(last_hit);
            busHit = gen(indLim2,1);
            bus(busHit,2) = 1; % Switch the bus to PQ node
            % deactivate the generator
            gen(indLim2,8) = 0;
            gen_b(indLim2) = 1;
            gen_a = ~gen_b;
            last_hit = indLim2;
        end
        
        if k > 0
            deltaLambda_loc = 0;
            bus(:,3) = P0;
            bus(:,4) = Q0;
        end
    end
end

if nbIter > maxIter
    success = 0;
end

% Pack the results
if System.dynData
    results = v2struct(success,nbItersPerStep,x,V,bus,deltaLambda_loc,gen_a,gen_b,last_hit,JacDyn,f_x,f_y,g_x,g_y);
else
    results = v2struct(success,nbItersPerStep,x,V,bus,deltaLambda_loc,gen_a,gen_b,last_hit,JacDyn);
end















