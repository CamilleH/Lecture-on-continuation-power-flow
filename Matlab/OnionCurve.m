% This file illustrates the onion curve in a two-bus system with one
% infinite bus and one load.
% In what follows, we study the closed-form solution of the load voltage
% magnitude at the load bus. We use normalized values i.e. V is actually
% V/E, Q is actually QX/E^2 and P is actually PX/E^2

close all;clc;clear;

% Two possible solutions for each pair (P,Q)
Vh = @(p,q) sqrt(1/2-q+sqrt(1/4-p.^2-q));
Vl = @(p,q) sqrt(1/2-q-sqrt(1/4-p.^2-q));
Plim = [0 0.8];
Qlim = [-0.4 0.3];

% List of values of p and q
nb_pts = 100;
P = linspace(Plim(1),Plim(2),nb_pts);
Q = linspace(Qlim(1),Qlim(2),nb_pts);
[PP,QQ] = meshgrid(P,Q);
% Compute the voltage solutions for all P and Q
all_Vh = zeros(nb_pts,nb_pts);
all_Vl = zeros(nb_pts,nb_pts);
for i = 1:nb_pts
    for j = 1:nb_pts
        all_Vh(i,j) = Vh(P(i),Q(j));
        all_Vl(i,j) = Vl(P(i),Q(j));
    end
end
% Transpose for plotting
all_Vh = all_Vh';
all_Vl = all_Vl';
% Remove non-physical solutions
all_Vh(imag(all_Vh)~=0) = NaN;
all_Vl(imag(all_Vl)~=0) = NaN;

% Cut part of the curves to show only load increases from zero load
all_Vh(PP<=-0.8/0.3*QQ) = NaN;
all_Vl(PP<=-0.8/0.3*QQ) = NaN;

% Plot the onion curve
onionfig = figure;
surf(PP,QQ,all_Vh);
hold on
surf(PP,QQ,all_Vl);
xlabel('$\frac{PX}{E^2}$','Interpreter','latex');
ylabel('$\frac{QX}{E^2}$','Interpreter','latex');
zlabel('V');
xlim(Plim);
ylim(Qlim);

%% Constant power factor loads
tanphi = [-0.3,0,0.4,1,10];
nb_tanphi = length(tanphi);

pvfig = figure;
nb_row = floor((nb_tanphi-1)/2)+1;
nb_col = 2;
for i = 1:nb_tanphi
    plotRectangle(tanphi(i),onionfig);
    Pvals = P;
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
    ylabel('V');
    title(sprintf('tanphi = %.1f',tanphi(i)));
end