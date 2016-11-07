% This file illustrates the effect of the power factor
close all;clc;clear;
% Defining the bounds of active and reactive power
Plim = [0 0.8];
Qlim = [-0.4 0.3];

% Two possible solutions for each pair (P,Q)
Vh = @(p,q) sqrt(1/2-q+sqrt(1/4-p.^2-q));
Vl = @(p,q) sqrt(1/2-q-sqrt(1/4-p.^2-q));

% Getting the onion figure
onionfig = get_onion_curve(Plim,Qlim,Vh,Vl);

% Constant power factor loads
tanphi = [-0.3,0,0.4,1,10];
nb_tanphi = length(tanphi);

pvfig = figure;
nb_row = floor((nb_tanphi-1)/2)+1;
nb_col = 2;
for i = 1:nb_tanphi
    plot_rectangle(tanphi(i),onionfig);
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
    figure(onionfig);
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