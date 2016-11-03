% This file shows how to determine operating points for different load
% models
clc;close all;clear;
% Defining the bounds of active and reactive power
Plim = [0 0.8];
Qlim = [-0.4 0.3];
Vlim = [0,1.4];

% Two possible solutions for each pair (P,Q)
Vh = @(p,q) sqrt(1/2-q+sqrt(1/4-p.^2-q));
Vl = @(p,q) sqrt(1/2-q-sqrt(1/4-p.^2-q));

%% Load model 1: constant power factor
% P = P0, Q = P*tanphi

% Getting the onion figure
onionfig = get_onion_curve(Plim,Qlim,Vh,Vl);

P0 = [0,0.2,0.4];
nb_p0 = length(P0);
tanphi = [-0.3,0,0.4,1,10];
nb_tanphi = length(tanphi);
colors = parula(nb_tanphi);
styles = ['-','--',':'];
all_lines = zeros(nb_tanphi,1);
all_lgd = cellstr(num2str(tanphi', 'tanphi=%.2f'));
for i = 1:nb_p0
    P_cur = P0(i);
    for j = 1:nb_tanphi
        Q_cur = P_cur*tanphi(j);
        x = P_cur*ones(2,1);
        y = Q_cur*ones(2,1);
        z = [Vlim(1);Vlim(2)];
        figure(onionfig);
        hold on
        % Create vertical bar corresponding to a certain active and
        % reactive power
        all_lines(j) = plot3(x,y,z,'Color',colors(j,:),'LineWidth',2);
    end
end
legend(all_lines,all_lgd);
title('Constant power factor load');

%% Load model 2: exponential load model
% P = z*P0*(V/V0)^alpha
% Q = z*Q0*(V/V0)^beta

% Getting the onion figure
onionfig2 = get_onion_curve(Plim,Qlim,Vh,Vl);

% Load model parameters
alpha = 1.5;
beta = 1.5;
P0 = 0.3;
Q0 = 0.2*P0;
V0 = 1;
P_exp = @(z,V) (z*P0*(V/V0).^alpha);
Q_exp = @(z,V) (z*Q0*(V/V0).^beta);

% Draw the load characteristics for different load demands z
nb_demand = 10;
demand = linspace(0,10,nb_demand);

% Range of load voltages
v_range = linspace(Vlim(1),Vlim(2),30);

all_load_charac = zeros(nb_demand,1);
demand_lgd = cellstr(num2str(demand', 'z=%.2f'));
for i = 1:nb_demand
    p_charac = P_exp(demand(i),v_range);
    q_charac = Q_exp(demand(i),v_range);
    figure(onionfig2);
    hold on;
    all_load_charac(i) = plot3(p_charac,q_charac,v_range);
end
legend(all_load_charac,demand_lgd);
title('Exponential load model');
