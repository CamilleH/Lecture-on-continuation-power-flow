function [F,G] = ch_calculateFGDyn(x,V,bus,gen_a,gen_b,Pg,Pwind,System)

indices = System.indices;
gen = System.gen;
% Update gen according to gen_b
gen(gen_b,8) = 0;
wind = System.wind;

if System.dynData
    gen_dyn = System.gen_dyn;
    ex = System.ex;
    param = System.param;
    Ybus_dyn = System.Ybus_dyn;
    Xdi = gen_dyn(:,6);
    Xdpi = gen_dyn(:,7);
    fs = param.fs;
    ws = fs*2*pi;
    Hi = gen_dyn(:,16);
    Mi = 2.*Hi./ws;
    Td0i = gen_dyn(:,9);
    Ka = ex(:,2);
    Tei = ex(:,5);
    Uref = ex(:,8);
    Ef_lim = ex(:,9);
else
    Ybus_stat = System.Ybus_stat;
end
%% Indices
% global indices
ng = indices.ng;
nbus = indices.nbus;

% indices of the Matpower arrays
GEN_BUS = 1;
PD = 3;
QD = 4;
PG = 2;

% Indices of the values in x
if System.dynData
    indX_delta = 1:ng;
    indX_omega = (ng+1):(2*ng);
    indX_Eqp = (2*ng+1):(3*ng);
    indX_Efp = (3*ng+1):(4*ng);
end

% Numerical indices of the sets a and b (AVRs that have not, and have,
% reached their limit, respectively);
gen_a_num = find(gen_a);
gen_b_num = find(gen_b);

%% Extract the values from x and y
if System.dynData
    delta = x(indX_delta);
    omega = x(indX_omega);
    Eqp = x(indX_Eqp);
    Ef = x(indX_Efp);
end

Va = angle(V);
Vm = abs(V);

%% Extracting the load and generation
SL = bus(:,PD) + 1i*bus(:,QD);
SWIND = zeros(nbus,1);
if ~isempty(wind)
    SWIND(wind(:,1)) = Pwind;
end
% Pg = gen(:,PG); % PG IS NOW A VARIABLE IN THE FUNCTION

if System.dynData
    % Voltages at the generator buses
    Vm_gen = Vm(gen(:,GEN_BUS));
    Va_gen = Va(gen(:,GEN_BUS));

    %% First calculate F as defined by the one-axis model plus AVR
    F1 = omega(2:ng);
    F2 = 1./Mi.*(Pg-Eqp.*Vm_gen./Xdpi.*sin(delta-Va_gen));
    F2(1) = []; % Taking away the row for the omega of the slack
    F3 = 1./Td0i.*(Ef - Xdi./Xdpi.*Eqp + (Xdi-Xdpi)./Xdpi.*Vm_gen.*cos(delta-Va_gen));
    
    % Here we calculate the equations for Ef, depending on whether the
    % exciter has hit its limits (for generators in gen_b) or not (in gen_a);
    Eq_a = 1./Tei.*(-Ef+Ka.*(Uref-Vm_gen));
    Eq_b = -Ef+Ef_lim;
    
    Eq_Ef_a =  Eq_a(gen_a_num);
    Eq_Ef_b = Eq_b(gen_b_num);
    
    F4 = Eq_Ef_a;
    
    F = [F1;F2;F3;F4];
    % F = [F2;F3;F4];
    
    %% Second we calculate G, all buses are PQ buses with Pgen = Qgen = 0
    EqpD = Eqp.*exp(1i*delta);
    V_dyn = [EqpD;V];
    Str = V_dyn.*conj(Ybus_dyn*V_dyn);
    Str(1:ng) = []; % no power balance at the generator buses
    G = [Eq_Ef_b;
        real(SL) + real(Str) - real(SWIND);
        imag(SL) + imag(Str) - imag(SWIND)];
else
    % Only the static equations are of interest (i.e. power equations)
    indGenBus = gen(:,1);
    on = gen_a;
    % generators off that have reached their Q limit
    offpq = gen_b ~= 0; 
    Str = V.*conj(Ybus_stat*V);
    F = [];
    DeltaP = real(SL) + real(Str) - real(SWIND);
    DeltaQ = imag(SL) + imag(Str) - imag(SWIND);
    % Subtracting active power generation from all generators that are
    % either on or off because they have reached their PQ limit
    DeltaP(indGenBus(on|offpq)) = DeltaP(indGenBus(on|offpq))-Pg(on|offpq);
    % Subtracting the reactive power generation of generators at the limit.
    DeltaQ(indGenBus(offpq)) = DeltaQ(indGenBus(offpq)) - gen(offpq,4); 
    pv = bus(:,2) == 2;
    pq = bus(:,2) == 1;
    % Adding the equations Vi = Vref for all sll generators, if any
    genSLL = gen_a==1 & gen_b==1;
    indBusSLL = gen(genSLL,1);
    if isempty(indBusSLL)
        eqVsll = [];
    else
        eqVsll = Vm(indBusSLL)-bus(indBusSLL,8);
    end
    G = [eqVsll;
        DeltaP(pv | pq);
        DeltaQ(pq)];
end
end
