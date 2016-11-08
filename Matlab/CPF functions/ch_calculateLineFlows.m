function [Sf,St] = ch_calculateLineFlows(V,System)

% global System

indices = System.indices;
Yt = System.Yt;
Yf = System.Yf;
Cf = System.Cf;
Ct = System.Ct;

ng = indices.ng;
nbus = indices.nbus;
ntot = ng+nbus;

% Extracting the voltages at the buses other than generator buses
Vstatic = V;

% Calculating the currents injected at both ends of the lines
If = Yf*Vstatic;
It = Yt*Vstatic;

% Getting the voltages at both ends
Vf = Cf*Vstatic;
Vt = Ct*Vstatic;

% Getting the injected powers
Sf = Vf.*conj(If);
St = Vt.*conj(It);