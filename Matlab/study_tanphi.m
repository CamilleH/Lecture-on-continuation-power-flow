% This file illustrates the effect of the power factor
close all;clc;clear;
% Defining the bounds of active and reactive power
Plim = [0 0.8];
Qlim = [-0.4 0.3];
Vlim = [0,1.4];

% Two possible solutions for each pair (P,Q)
Vh = @(p,q) sqrt(1/2-q+sqrt(1/4-p.^2-q));
Vl = @(p,q) sqrt(1/2-q-sqrt(1/4-p.^2-q));

%% Load model: constant power factor
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

%% All solutions for constant power factor
% Getting the onion figure
onionfig2 = get_onion_curve(Plim,Qlim,Vh,Vl);

% Constant power factor loads
tanphi = [-0.3,0,0.4,1,10];
nb_tanphi = length(tanphi);

pvfig = figure;
nb_row = floor((nb_tanphi-1)/2)+1;
nb_col = 2;
for i = 1:nb_tanphi
    plot_rectangle(tanphi(i),onionfig2);
    Pvals = linspace(Plim(1),Plim(2),100);
    Qvals = Pvals*tanphi(i);
    Vl_sol = Vl(Pvals,Qvals);
    Vh_sol = Vh(Pvals,Qvals);
    % Remove complex solutions
    Vl_sol(imag(Vl_sol)~=0) = NaN;
    Vh_sol(imag(Vh_sol)~=0) = NaN;
    % Bounding Q
    Qvals(Qvals<Qlim(1)|Qvals>Qlim(2)) = NaN;
    % Plotting the curves on the onion curve
    figure(onionfig2);
    hold on;
    plot3(Pvals,Qvals,Vl_sol,'r','LineWidth',2);
    plot3(Pvals,Qvals,Vh_sol,'r','LineWidth',2);
    % Plotting the curves in the PV plane
    figure(pvfig);
    subplot(nb_row,nb_col,i);
    plot(Pvals,Vl_sol,'k');
    hold on;
    plot(Pvals,Vh_sol,'k');
    xlabel('$\frac{PX}{E^2}$','Interpreter','latex');
    ylabel('$\frac{V}{E}$','Interpreter','latex');
    title(sprintf('tanphi = %.1f',tanphi(i)));
end