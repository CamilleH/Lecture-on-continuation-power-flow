% This file compares PV and QV curves
close all;clc;clear;
% Defining the bounds of active and reactive power
Plim = [0 0.8];
Qlim = [-0.4 0.3];
Vlim = [0 1.4];

% Two possible solutions for each pair (P,Q)
Vh = @(p,q) sqrt(1/2-q+sqrt(1/4-p.^2-q));
Vl = @(p,q) sqrt(1/2-q-sqrt(1/4-p.^2-q));

% Getting the onion figure
onionfig = get_onion_curve(Plim,Qlim,Vh,Vl);

% Choosing one initial operating point
Pop = 0.2;
Qop = 0.1;

% Draw the initial point
x = Pop*ones(2,1);
y = Qop*ones(2,1);
z = [Vlim(1);Vlim(2)];
figure(onionfig);
hold on
plot3(x,y,z,'r','LineWidth',2);
%% PV curve
% To draw the PV curve, we must assume a certain relationship between P and
% Q, for example, constant power factor
tanphi = 0.5;
% We draw the plane corresponding to this tanphi
plot_rectangle(tanphi,onionfig);
% We draw the set of operating points corresponding to this tangent phi
P_pv = linspace(Pop,Plim(2),100);
Q_pv = tanphi*P_pv;
V_lowpv = Vl(P_pv,Q_pv);
V_highpv = Vh(P_pv,Q_pv);
V_lowpv(imag(V_lowpv) ~= 0) = NaN;
V_highpv(imag(V_highpv) ~= 0) = NaN;
figure(onionfig);
hold on;
plot3(P_pv,Q_pv,V_lowpv,'k','LineWidth',2);
plot3(P_pv,Q_pv,V_highpv,'k','LineWidth',2);
% Projection to the PV plane
pvfig = figure;
hold on;
plot(P_pv,V_lowpv,'k','LineWidth',2);
plot(P_pv,V_highpv,'k','LineWidth',2);
xlabel('P');
ylabel('V');

%% QV curve
plot_rectangle_QV(Pop,onionfig);
% Plot the set of operating points corresponding to the QV curve
Q_qv = linspace(-0.08,0.3,100);
P_qv = Pop*ones(1,length(Q_qv));
V_lowqv = Vl(P_qv,Q_qv);
V_highqv = Vh(P_qv,Q_qv);
V_lowqv(imag(V_lowqv)~=0) = NaN;
V_highqv(imag(V_highqv)~=0) = NaN;
% Find the reactive power injection needed to get a voltage of 0.95
idx095 = find(V_highqv > 0.945 & V_highqv < 0.955);
Q095 = Q_qv(idx095(1));
% Draw on the onion curve
figure(onionfig);
hold on
plot3(P_qv,Q_qv,V_highqv,'b','LineWidth',2);
plot3(P_qv,Q_qv,V_lowqv,'b','LineWidth',2);
% Draw the QV curve
qvcurve = figure;
hold on;
plot(V_highqv,-Q_qv,'k','LineWidth',2);
plot(V_lowqv,-Q_qv,'k','LineWidth',2);
curop = plot([min(V_lowqv);max(V_highqv)],-Qop*ones(2,1),'--','LineWidth',2);
op095 = plot([min(V_lowqv);max(V_highqv)],-Q095*ones(2,1),'--','LineWidth',2);
xlabel('V');
ylabel('Reactive power injection');
legend([curop,op095],{'Current Q injection','Q injection needed for V=0.95'});
