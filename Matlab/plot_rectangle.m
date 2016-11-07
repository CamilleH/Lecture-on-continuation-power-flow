function plot_rectangle(tanphi,hfig)
Plim = [0,0.8];
Qlim = [-0.4,0.3];
Vlim = [0,1.4];
rectP = [Plim(1) Plim(1);Plim(2) Plim(2)];
rectQ = tanphi*rectP;
rectV = [Vlim(1) Vlim(2);Vlim(1) Vlim(2)];
% Bounding Q
idxQ_toolow = rectQ<Qlim(1);
idxQ_toohigh = rectQ>Qlim(2);
rectP(idxQ_toolow) = Qlim(1)/tanphi;
rectQ(idxQ_toolow) = Qlim(1);
rectP(idxQ_toohigh) = Qlim(2)/tanphi;
rectQ(idxQ_toohigh) = Qlim(2);
% Draw rectangle
recColor = zeros(2,2,3);
recColor(:,:,:) = 0.5;
figure(hfig);
hold on
srec1 = surf(rectP,rectQ,rectV,recColor);
alpha(srec1,0.3);
end