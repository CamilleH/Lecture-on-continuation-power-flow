function [xpred,Vpred,lambdapred,t] = ch_predictDyn(System,x,V,lambda,gen_a,gen_b,dirP,dirQ,step,bus,k)

indices = System.indices;

if nargin == 9
    k = 0;
end

%% Indices
% global indices
ng = indices.ng;
nbus = indices.nbus;

gen_a_num = find(gen_a);
gen_b_num = find(gen_b);
ngen_a = length(gen_a_num);

if System.dynData
    n_x = length(x)-2-length(gen_b_num);
    n_y = 2*nbus+length(gen_b_num);
    n_z = n_x+n_y;
    indZ_U = 4*ng-2+nbus+(1:nbus);
    indF_P = 4*ng-2+(1:nbus);
    indF_Q = 4*ng-2+nbus+(1:nbus);
    indX = 1:(4*ng-2);
    indY = (4*ng-1):(4*ng+2*nbus-2);
    indVa = true*ones(nbus,1);
    indVm = true*ones(nbus,1);
else
    n_y = 2*nbus;
    n_z = n_y;
    pq = bus(:,2) == 1;
    pv = bus(:,2) == 2;
    slack = bus(:,2) == 3;
    indVa = ~slack;
    indVm = ~(slack | pv);
    indVm_num = find(indVm);
    indZ_U = nbus+1:2*nbus;
    eqPremove = ~(pv | pq);
    eqQremove = ~pq;
    indF_P = 1:nbus;
    indF_Q = (nbus+1):2*nbus;
    indX = [];
    indY = 1:(sum(indVa)+sum(indVm));
end

%% Building the system of equations

if size(dirP,2) > 1
    dirP = dirP.';
    dirQ = dirQ.';
end

[JacPF,~] = ch_calculateJacDyn(x,V,gen_a,gen_b,System,bus);
F_z = JacPF;
F_lambda = zeros(n_z,1);
F_lambda(indF_P) = dirP;
F_lambda(indF_Q) = dirQ;

e = zeros(1,n_z+1);
if k == 0
    e(end) = 1;
elseif k == -1
    e(end) = -1;
else
    indkine = indZ_U(k);
    e(indkine) = -1;
end

if ~System.dynData
    % order important here since the first line removes elements at the end
    % of the matrix, and the second at the beginning
    F_lambda(nbus+find(eqQremove)) = [];
    F_lambda(eqPremove) = [];
    e(nbus+find(eqQremove)) = [];
    e(eqPremove) = [];
end

A = [F_z F_lambda;
    e];
B = [zeros(size(F_z,1),1);
    1];

deltaXL = A\B;
t=deltaXL;
t=t/norm(t);
deltaXL = step*deltaXL;
deltaX = deltaXL(indX);
deltaY = deltaXL(indY);
deltaLambda = deltaXL(end);

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
Va = angle(V);
Vm = abs(V);
Va(indVa) = Va(indVa) + deltaY(indVa);
Vm(indVm) = Vm(indVm) + deltaY(sum(indVa)+(1:sum(indVm)));

Vpred = Vm.*exp(1i*Va);
xpred = x;
lambdapred = lambda+deltaLambda;
end