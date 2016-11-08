function [k,slopeVmin] = ch_choosePara(System,Settings_ch,x,V,dirP,dirQ,gen_a,gen_b,bus,varargin)

cpfSettings = Settings_ch.cpfSettings;
slopeLim = cpfSettings.slopeLim1;

flagForce = length(varargin);

%% Indices
ng = length(gen_a);%indices.ng;
nbus = size(bus,1);%indices.nbus;

gen_a_num = find(gen_a);
gen_b_num = find(gen_b);
ngen_a = length(gen_a_num);

if System.dynData
    n_x = length(x)-2-length(gen_b_num);
    n_y = 2*nbus+length(gen_b_num);
    n_z = n_x+n_y;
    indF_P = 4*ng-2+length(gen_b_num)+(1:nbus);
    indF_Q = 4*ng-2+nbus+(1:nbus);
    indZ_Eqp = (2*ng+1):(3*ng);
    indZ_Eqp = indZ_Eqp-2;
    % indY_U = length(gen_b_num)+nbus+(1:nbus);
    indZ_Vm = 4*ng-2+nbus+(1:nbus);
else
    n_x = 0;
    n_y = 2*nbus;
    pq = bus(:,2) == 1;
    pv = bus(:,2) == 2;
    slack = bus(:,2) == 3;
    eqPremove = ~(pv | pq);
    eqQremove = ~pq;
    indVa = ~slack;
    indVm = ~(pv | slack);
    indVm_num = find(indVm);
    indF_P = 1:nbus;
    indF_Q = nbus+1:2*nbus;
    indZ_Vm = sum(indVa)+(1:sum(indVm));
    indZ_Eqp = [];
end
%% Calculate the tangent
if size(dirP,2) > 1
    dirP = dirP.';
    dirQ = dirQ.';
end

[JacPF,~] = ch_calculateJacDyn(x,V,gen_a,gen_b,System,bus);
F_z = JacPF;
F_lambda = zeros(n_x+n_y,1);
F_lambda(indF_P) = dirP;
F_lambda(indF_Q) = dirQ;
if ~System.dynData
    % order important here since the first line removes elements at the end
    % of the matrix, and the second at the beginning
    F_lambda(nbus+find(eqQremove)) = [];
    F_lambda(eqPremove) = [];
end

if size(F_z,1) ~= size(F_lambda,1)
    keyboard;
end

dz_dlambda = -F_z\F_lambda;

dEqp_dlambda = dz_dlambda(indZ_Eqp);
dU_dlambda = dz_dlambda(indZ_Vm);

if size(dEqp_dlambda,2) > 1
    dEqp_dlambda = dEqp_dlambda.';
    dU_dlambda = dU_dlambda.';
end
% dV_dlambda = [dEqp_dlambda;dU_dlambda];
dV_dlambda = dU_dlambda;
[slopeVmin,indMin] = min(dV_dlambda);

if (abs(slopeVmin) > slopeLim) || flagForce
    % we choose the critical voltage either when its derivative is large
    % enough, or when we want to force using another parameter than lambda
    if System.dynData
        k = indMin;
    else
        % IN this case, we look only at a subset of all voltages, since the
        % voltages at the slack and pv buses are kept constant
        k = indVm_num(indMin);
    end
%     keyboard;
else
    k = 0;
end

end