function plot_rectangle_QV(Pop,hfig)
% This function draws a plane going through Pop, and corresponding
% to changes in Q only.
Qlim = [-0.4,0.3];
Vlim = [0,1.4];
rectP = Pop*ones(2);
rectQ = [Qlim(1) Qlim(1);Qlim(2) Qlim(2)];
rectV = [Vlim(1) Vlim(2);Vlim(1) Vlim(2)];
% Draw rectangle
recColor = zeros(2,2,3);
recColor(:,:,:) = 0.5;
figure(hfig);
hold on
srec1 = surf(rectP,rectQ,rectV,recColor);
alpha(srec1,0.5);
end