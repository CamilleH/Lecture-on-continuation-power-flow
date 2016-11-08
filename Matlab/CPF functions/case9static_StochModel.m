function stochMod = case9static_StochModel(caseName,wpNormal)
% Describes the stochastic model used fo the system

% Loading the parameters of the wind farm (installed cap)
mpc = openCase('case9static');
if isempty(mpc.wind)
    WPcap = [];
else
    WPcap = mpc.wind(:,4); % installed wind power capacity
end

% Loading the settings
Settings_ch = ch_settings_init();

% Most likely point
as=25.'; % a parameter for the beta distribution
bs=25; % b parameter for the beta distribution
nwp = length(WPcap);
sigWP = eye(nwp);

if strcmp(caseName,'load5')
    muLoad = mpc.bus(5,3);
elseif strcmp(caseName,'load6')
    muLoad = mpc.bus(6,3);
elseif strcmp(caseName,'load8')
    muLoad = mpc.bus(8,3);
elseif strcmp(caseName,'loads579')
    muLoad = mpc.bus([5 7 9],3);
end
nload = length(muLoad);
SigLoad = 0.5^2*eye(nload); % standard deviation of the load

nBeta = length(WPcap);
nNorm = length(muLoad);

withU = Settings_ch.bifSettings.inclContGen;
if withU
    % Getting optimal direction for u
    [du_star_dzeta,u_opt] = getOptimDudzeta('case9static',caseName);
    u_opt = 1;
    Gamma=du_star_dzeta;
    sigBeta=sqrt(WPcap.^2.*as.*bs./((as+bs).^2.*(as+bs+1)));
    Sig_Theta=[diag(sigBeta.^2) zeros(nBeta,nNorm);
        zeros(nNorm,nBeta) SigLoad];
    Sig_u=diag(diag(Gamma*Sig_Theta*Gamma'));
    indU = 1:length(u_opt);
else
    indU = [];
    u_opt = [];
    Sig_u = [];
end

% indBeta = 1; % indices of WP in lambda, beta distributed
% indNorm = 2; % indices of the lambda in zeta, normal distributed

stochMod = ModelUNormZetaBetaNorm(indU,WPcap,as,bs,...
    muLoad,SigLoad,u_opt,Sig_u,sigWP,wpNormal);