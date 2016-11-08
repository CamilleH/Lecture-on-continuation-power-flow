function deltaQ = computeGenQ(System,V,bus)
% This function computes the reactive power production from the generators
% necessary to fulfill the power flow equations
Ybus_stat = System.Ybus_stat;
Str = V.*conj(Ybus_stat*V);
SL = bus(:,3) + 1i*bus(:,4);
deltaQ = imag(SL) + imag(Str);
end
