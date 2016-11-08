classdef (Abstract) ModelUncert < handle
% Abstact class used as an interface with the different uncertainty models

    methods(Abstract = true)
        f(obj,lambda); % importance function
        Df(obj,lambda); % Derivative of the importance function
    end
    
end