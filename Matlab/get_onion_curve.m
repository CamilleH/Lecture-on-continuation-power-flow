function onionfig = get_onion_curve(Plim,Qlim,Vh,Vl)
% This file illustrates the onion curve in a two-bus system with one
% infinite bus and one load.
% In what follows, we study the closed-form solution of the load voltage
% magnitude at the load bus. We use normalized values i.e. V is actually
% V/E, Q is actually QX/E^2 and P is actually PX/E^2

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
zlabel('$\frac{V}{E}$','Interpreter','latex');
xlim(Plim);
ylim(Qlim);
view(39,50);

% Add a thick black line for the maximum loadability points
Vmax = @(phi) (1./(sqrt(2)*sqrt(1+sin(phi))));
Pmax = @(phi) (cos(phi)./(2*(1+sin(phi))));
Qmax = @(phi) (sin(phi)./(2*(1+sin(phi))));
tanphi = linspace(-0.3/0.8,10,30);
phi = atan(tanphi);
figure(onionfig);
hold on;
% Plot the maximum loadability limits
plot3(Pmax(phi),Qmax(phi),Vmax(phi),'b','LineWidth',2);
% Projection on the (P,Q) plane
plot3(Pmax(phi),Qmax(phi),zeros(1,length(phi)),'b:','LineWidth',2);
