function [JacPF,JacDyn,f_x,f_y,g_x,g_y] = ch_calculateJacDyn(x,V,gen_a,gen_b,System,bus,flagNoJacDyn)
% Calculate the dynamic Jacobian according to As = fx-fy*(gy)^-1*gx
%    with x = [delta...,omega...,Eqp...,Efp(gen_a)...,Efp(gen_b),...]
%    and  y = [theta...,U...]
%
%    The first gen is assumed to be the slack

if nargin == 6
    flagNoJacDyn = 0;
end

indices = System.indices;
gen = System.gen;

if System.dynData
    param = System.param;
    gen_dyn = System.gen_dyn;
    ex = System.ex;
else
    flagNoJacDyn = 1;
end

%% Ybus and parameters
% Ybus_stat = System.Ybus_stat;
if System.dynData
    Ybus_dyn = System.Ybus_dyn;
    Xdi = gen_dyn(:,6);
    Xdpi = gen_dyn(:,7);
    Di = gen_dyn(:,17);
    fs = param.fs;
    ws = fs*2*pi;
    Hi = gen_dyn(:,16);
    Mi = 2.*Hi./ws;
    Td0i = gen_dyn(:,9);
    Ka = ex(:,2);
    Tei = ex(:,5);
    % Uref = ex(:,8);
    % Ef_lim = ex(:,9);
else
    Ybus_stat = System.Ybus_stat;
end

%% Indices

% global indices
ng = indices.ng;
nbus = indices.nbus;

% indices of the Matpower arrays
GEN_BUS = 1;
% PD = 3;
% QD = 4;

% Indices of the values in x and y
if System.dynData
    indX_delta = 1:ng;
    indX_omega = (ng+1):(2*ng);
    indX_Eqp = (2*ng+1):(3*ng);
    indX_Efp = (3*ng+1):(4*ng);
end

% Indices of the buses at which the generators are connected (not Eq but U)
ind_gen = gen(:,GEN_BUS).';

% Numerical indices of the sets a and b (AVRs that have not, and have,
% reached their limit, respectively);
gen_a_num = find(gen_a);
gen_b_num = find(gen_b);
ngen_b = length(gen_b_num);

% Find a generator that is both in aset and bset (for SLL). The sets below
% do not have the generator in common (if any)
gen_b_NOsll = ~gen_a;
gen_a_NOsll = ~gen_b;
gen_b_NOsll_num = find(gen_b_NOsll);
gen_a_NOsll_num = find(gen_a_NOsll);

% indices to remove
if System.dynData
    indF_remove = [1;ng+1;3*ng+gen_b_NOsll_num];
    indG_remove = gen_a_NOsll_num;
    indX_remove = [1;ng+1;3*ng+gen_b_NOsll_num];
    indY_remove = gen_a_num; % we remove the possible generator in common from the Yis
end

%% Extract the values from x and y
if System.dynData
    delta = x(indX_delta);
    % omega = x(indX_omega);
    Eqp = x(indX_Eqp);
    % Ef = x(indX_Efp);
end

Va = angle(V);
Vm = abs(V);

% Voltages at the generator buses
Vm_gen = Vm(ind_gen);
Va_gen = Va(ind_gen);

if System.dynData
    %% fx fy
    %  fx = [f_delta_x f_delta_y;
    %        f_omega_x f_omega_y;
    %        f_Eq_x    f_Eq_y;
    %        F_Ef_x    f_Ef_y]
    
    f_delta_x = sparse(1:ng,indX_omega,1,ng,4*ng);
    
    % f_omega_x = [f_om_delta f_om_om f_om_Eq f_om_Ef];
    val_om_delta = -Eqp.*Vm_gen.*cos(delta-Va_gen)./(Mi.*Xdpi);
    val_om_om = -Di./Mi;
    val_om_Eq = -Vm_gen.*sin(delta-Va_gen)./(Mi.*Xdpi);
    f_omega_x = sparse([1:ng,1:ng,1:ng],[indX_delta,indX_omega,indX_Eqp],[val_om_delta;val_om_om;val_om_Eq],ng,4*ng);
    % f_omega_x = sparse([1:ng,1:ng],[indX_delta,indX_Eqp],[val_om_delta;val_om_Eq],ng,3*ng);
    
    % f_Eq_x
    val_Eq_delta = -(Xdi-Xdpi).*Vm_gen.*sin(delta-Va_gen)./(Td0i.*Xdpi);
    val_Eq_Eq = -Xdi./(Td0i.*Xdpi);
    val_Eq_Ef = 1./(Td0i);
    f_Eq_x = sparse([1:ng,1:ng,1:ng],[indX_delta,indX_Eqp,indX_Efp],[val_Eq_delta;val_Eq_Eq;val_Eq_Ef],ng,4*ng);
    
    % f_Ef_x
    val_Ef_Ef = -1./Tei;
    f_Ef_x = sparse(1:ng,indX_Efp,val_Ef_Ef,ng,4*ng);
    % f_Ef_x(gen_b_sll_num,:) = []; % Taking away the exciters that have reached their limits
    
    % f_x
    f_x = [f_delta_x;f_omega_x;f_Eq_x;f_Ef_x];
    % f_x = [f_omega_x;f_Eq_x;f_Ef_x];
    % f_x(:,3*ng+gen_b_num) = []; % Taking away the Ef that have reached their limits from the state variables
    % f_x(ng+1,:) = []; % Taking away the omega equation for the slack
    % f_x(:,ng+1) = []; % Taking away the omega for the slack from the state variables
    % f_x(1,:) = []; % Taking away the row for the delta of the slack
    % f_x(:,1) = []; % Taking away the column for the delta of the slack
    
    % f_y: don't forget the Ef that have become algebraic variables because the
    % corresponding AVRs have reached their limit.
    % f_delta_y = 0
    % f_om_y, the first ng columns/rows correspond to Ef.
    val_om_Va = Eqp.*Vm_gen./(Mi.*Xdpi).*cos(delta-Va_gen);
    val_om_Vm = -Eqp./(Mi.*Xdpi).*sin(delta-Va_gen);
    f_om_y = sparse([1:ng,1:ng],[ng+ind_gen,ng+nbus+ind_gen],[val_om_Va;val_om_Vm],ng,2*nbus+ng);
    
    % f_Eq_y
    val_Eq_Ef = 1./(Td0i); % case where some AVRs have hit their limit
    val_Eq_Va = (Xdi-Xdpi)./(Td0i.*Xdpi).*Vm_gen.*sin(delta-Va_gen);
    val_Eq_Vm = (Xdi-Xdpi)./(Td0i.*Xdpi).*cos(delta-Va_gen);
    f_Eq_y = sparse([1:ng,1:ng,1:ng],[1:ng,ng+ind_gen,ng+nbus+ind_gen],[val_Eq_Ef;val_Eq_Va;val_Eq_Vm],ng,2*nbus+ng);
    
    % f_Ef_y
    val_Ef_Ef = 1./Tei;
    val_Ef_Vm = -Ka./Tei;
    f_Ef_y = sparse([1:ng,1:ng],[1:ng,ng+nbus+ind_gen],[val_Ef_Ef;val_Ef_Vm],ng,2*nbus+ng);
    % f_Ef_y(gen_b_sll_num,:) = []; % Taking away the exciters that have reached their limits
    
    % f_y
    f_y = [sparse(zeros(ng,2*nbus+ng));f_om_y;f_Eq_y;f_Ef_y];
    % f_y = [f_om_y;f_Eq_y;f_Ef_y];
    % f_y(ng+1,:) = []; % Taking away the omega equation for the slack
    % f_y(1,:) = []; % taking away the first row that corresponds to the angle of the slack bus
    % f_y(:,gen_a_num) = []; % Taking away the columns corresponding to the AVR that have not reached their limits
    
    %% gx gy
    EqpD = Eqp.*exp(1i*delta);
    V_dyn = [EqpD;V];
    % The first rows and columns in Ybus_dyn are the generators' internal nodes
    % The first voltages are these at the generators' internal nodes.
    [dSbus_dVm, dSbus_dVa] = dSbus_dV(Ybus_dyn, V_dyn);
    % [dSbus_dVm2, dSbus_dVa2] = ch_calcJacTens(Ybus_dyn, V_dyn);
    
    % gx = [dEflim/dDelta dEflim/dOmega dEflim/dEqp dEflim/dEf;
    %       dP/dDelta dP/dOmega dP/dEqp dP/dEf;
    %       dQ/dDelta dQ/dOmega dQ/dEqp dQ/dEf];
    % g_Ef = zeros(ng,4*ng)); % None of the Ef that are included in g are in x (all in y)
    g_Ef = sparse(1:ng,indX_Efp,-ones(ng,1),ng,4*ng);
    % g_Ef(gen_a_num,:) = []; % taking away the exciters that have not reached their limits
    g_x = [g_Ef;
        real(dSbus_dVa((ng+1):(ng+nbus),1:ng)),sparse(zeros(nbus,ng)),real(dSbus_dVm((ng+1):(ng+nbus),1:ng)),sparse(zeros(nbus,ng));
        imag(dSbus_dVa((ng+1):(ng+nbus),1:ng)),sparse(zeros(nbus,ng)),imag(dSbus_dVm((ng+1):(ng+nbus),1:ng)),sparse(zeros(nbus,ng))];
    % g_x = [g_Ef;
    %        real(dSbus_dVa((ng+1):(ng+nbus),1:ng)),real(dSbus_dVm((ng+1):(ng+nbus),1:ng)),sparse(zeros(nbus,ng));
    %        imag(dSbus_dVa((ng+1):(ng+nbus),1:ng)),imag(dSbus_dVm((ng+1):(ng+nbus),1:ng)),sparse(zeros(nbus,ng))];
    % g_x(:,3*ng+gen_b_num) = []; % Taking away the Ef of the AVRs that have reached their limit
    % g_x(:,ng+1) = []; % Taking away the omega equation for the slack
    % g_x(:,1) = []; % Taking away the angle of the slack bus
    
    % gy = [dEflim/dEf dEflim/dVa dEflim/dVm;
    %       dP/dEf dP/dVa dP/dVm;
    %       dQ/dEf dQ/dVa dQ/dVm];
    g_Eflim_y = sparse(gen_b_num,gen_b_num,-ones(ngen_b,1),ng,ng+2*nbus);
    % g_P_Ef = real(dSbus_dVm((ng+1):(ng+nbus),1:ng));
    % g_Q_Ef = imag(dSbus_dVm((ng+1):(ng+nbus),1:ng));
    % g_P_Ef(:,gen_a_num) = []; % Taking away the columns of the Ef that are still state variables in x
    % g_Q_Ef(:,gen_a_num) = [];
    g_P_Ef = zeros(nbus,ng);
    g_Q_Ef = zeros(nbus,ng);
    g_y = [g_Eflim_y;
        g_P_Ef real(dSbus_dVa((ng+1):(ng+nbus),(ng+1):(ng+nbus))) real(dSbus_dVm((ng+1):(ng+nbus),(ng+1):(ng+nbus)));
        g_Q_Ef imag(dSbus_dVa((ng+1):(ng+nbus),(ng+1):(ng+nbus))) imag(dSbus_dVm((ng+1):(ng+nbus),(ng+1):(ng+nbus)))];
    
    %% removing indices
    f_x(indF_remove,:) = [];
    f_y(indF_remove,:) = [];
    g_x(indG_remove,:) = [];
    g_y(indG_remove,:) = [];
    f_x(:,indX_remove) = [];
    g_x(:,indX_remove) = [];
    f_y(:,indY_remove) = [];
    g_y(:,indY_remove) = [];
    
    
    %% Jacobians
    if ~flagNoJacDyn
        temp = g_y\g_x;
        JacDyn = f_x - f_y*temp;
    else
        JacDyn = [];
    end
    JacPF = [f_x f_y;
        g_x g_y];
else
    % generator types
    pv = bus(:,2) == 2;
    pq = bus(:,2) == 1;
    % if sll
    indGenSLL = gen_a == 1 & gen_b == 1;
    % Just static data (power flow equations only)
    [dSbus_dVm, dSbus_dVa] = dSbus_dV(Ybus_stat, V);
    dP_dVa = real(dSbus_dVa);
    dP_dVm = real(dSbus_dVm);
    dQ_dVa = imag(dSbus_dVa);
    dQ_dVm = imag(dSbus_dVm);
    JacPF = [dP_dVa(pv | pq,pv | pq) dP_dVm(pv | pq,pq);dQ_dVa(pq,pv | pq) dQ_dVm(pq,pq)];
    % Add equations of V = Vref for sll generators
    indBusSLL = System.gen(indGenSLL,1);
    % check whether the corresponding bus is a PQ bus (as it should)
    if nnz(pq(indBusSLL) == 0) > 0
        keyboard;
    end
    % numerical index of the SLL Vm among all pq indices
    indBusSLLinPQ = find(ismember(find(pq),indBusSLL)); 
    indJacVmSLL = nnz(pv)+nnz(pq)+indBusSLLinPQ;
    vrefsll = zeros(length(indBusSLL),nnz(pq)+nnz(pv)+nnz(pq));
    vrefsll(:,indJacVmSLL) = eye(length(indBusSLL));
    JacDyn = [vrefsll;JacPF];
    JacPF = JacDyn;
    f_x = [];
    f_y = [];
    g_x = [];
    g_y = [];
end
end
