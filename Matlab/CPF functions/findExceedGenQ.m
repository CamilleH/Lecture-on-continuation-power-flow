function indLim = findExceedGenQ(deltaQ,gen,genToCheck,varargin)
% This function returns the indices of all generators in the set genToCheck
% that exceed their limits
if isempty(genToCheck)
    keyboard;
end
if isempty(varargin)
    mode = 'below';
else
    mode = varargin{1};
end
Qg = deltaQ(gen(genToCheck,1));
Qmax = gen(genToCheck,4);
if strcmp(mode,'below')
    tol = -1e-8;
elseif strcmp(mode,'exceed')
    tol = 1e-8;
else
    keyboard;
end
diffQ = (Qg - Qmax) > tol;
indLim = genToCheck(diffQ);