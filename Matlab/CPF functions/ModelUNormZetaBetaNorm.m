classdef ModelUNormZetaBetaNorm < ModelUncert
    properties
        indU; % indices of u in lambda
        indBeta; % indices of the beta distributed variables in lambda
        indNorm; % indices of the Gaussian distributed variables in lambda
        WPcap; % normalization factor for the beta (ex: Wcap)
        nwind;
        nload;
        nu;
%         Cp; % Capacity factors of the wind power
        % Coefficient by which the forecast errors should be multiplied in 
        % order for the prod to always be negative and smaller than the 
        % installed cap
%         Kfore; 
        % Parameters of the beta distribution for the forecast errors
        a;
        b;
        muBeta;
        sigBeta;
        % Parameters of the normal distribution (could be
        % multidimensionals)
        muNorm;
        sigNorm;
        % parameters for u
        muU;
        sigU;
        % Covariance matrix of all zeta (beta first and then normal)
        muAll;
        muX;
        sigX;
        muXwp;
        muWl;
        zetaX;
        % Possibly: covariance matrix for the wind power
        sigWP;
        wpNormal;
        modeStudentt = 0;
        muWPis;
        % For decorrelated variables
        Mwp;
        Ml;
        Mchol; % such that Sigma = M M^T, this is for Sigma in zeta space
        % Flag to consider u or not
        flagU = 1;
        flagDecorr = 0;
        flagQuadr = 1; % flag for using (x-\mu)\Sig^{-1}(x-\
%         sigAll;
    end
    
    methods
        function obj = ModelUNormZetaBetaNorm(indU,WPcap,a,b,mu,sig,muU,sigU,sigWP,wpNormal)
            obj.indU = indU;
%             obj.indBeta = indBeta;
%             obj.indNorm = indNorm;
            obj.WPcap = WPcap;
            obj.nwind = length(WPcap);
            obj.nload = length(mu);
            obj.nu = length(obj.indU);
            obj.indBeta = obj.nu+1:obj.nu+obj.nwind;
            obj.indNorm = obj.nu+obj.nwind+1:obj.nu+obj.nwind+obj.nload;
%             obj.Cp = Cp;
%             obj.Kfore = Kfore;
            obj.a = a;
            obj.b = b;
            obj.muNorm = mu;
            obj.sigNorm = sig; % Covariance matrix
            obj.muU = muU;
            obj.sigU = sigU; % Covariance matrix
            obj.sigWP = sigWP;
            obj.computeMuSig();
            obj.muWPis = zeros(size(obj.sigWP,1),1);
            obj.Mwp = chol(sigWP).';
            obj.Ml = chol(obj.sigNorm).';
            obj.Mchol = [obj.Mwp zeros(obj.nwind,obj.nload);
                zeros(obj.nload,obj.nwind) obj.Ml];
            obj.wpNormal = wpNormal;
        end
        
        function updateToNewStochModel(obj)
            obj.nwind = length(obj.WPcap);
            obj.nload = length(obj.muNorm);
            obj.nu = length(obj.indU);
            obj.computeMuSig();
            obj.Mwp = chol(obj.sigWP).';
            obj.Ml = chol(obj.sigNorm).';
            obj.Mchol = [obj.Mwp zeros(obj.nwind,obj.nload);
                zeros(obj.nload,obj.nwind) obj.Ml];
        end
        
        function stMod_copy = copy(obj)
             stMod_copy = ModelUNormZetaBetaNorm(obj.indU,obj.WPcap,obj.a,obj.b,...
                 obj.muNorm,obj.sigNorm,obj.muU,obj.sigU,obj.sigWP,obj.wpNormal);
        end
        
        function computeMuSig(obj)
            % Computations of other parameters
            obj.muBeta = obj.WPcap.*obj.a./(obj.a+obj.b);
            obj.muAll = [obj.muBeta;obj.muNorm];
            obj.sigBeta = sqrt(obj.WPcap.*obj.a.*obj.b./((obj.a+obj.b).^2.*(obj.a+obj.b+1))); % standard deviation
            % In the X = tilde(zeta) space
            obj.muXwp = zeros(obj.nwind,1);
            obj.muX = [obj.muXwp;obj.muNorm];
            obj.sigX = [obj.sigWP zeros(obj.nwind,obj.nload);
                zeros(obj.nload,obj.nwind) obj.sigNorm];
        end
        
        function changeLoads(obj,mu,sig)
            obj.muNorm = mu;
            obj.sigNorm = sig;
            obj.computeMuSig();
        end
        
        function changeWP(obj,mu,sig)
            obj.muWPis = mu;
            if nargin == 3
                obj.sigWP = sig;
            end
            obj.computeMuSig();
        end
        
        %% Functions to link forecast errors and wind power productions
        function x = forecasttoprod(obj,z)
            x = obj.WPcap.*z;
        end
        
        function [z,varargout] = prodtoforecast(obj,x)
            z = x./obj.WPcap;
            if nargout > 1
                % then we also compute the derivative of z wrt x
                dz_dx = 1./obj.WPcap;
                varargout{1} = dz_dx;
            end
        end
        
        
        %% Values of the uncertainty for a given zeta
        function value = fnorm(obj,x)
            if size(x,2) > 1
                x = x.';
            end
            
            if ~obj.modeStudentt
                value = mvnpdf(x,obj.muNorm,obj.sigNorm);
            else
                % We use the Student's t distribution instead, for IS
                value = mvtpdf(obj.sigNorm\(x-obj.muNorm),eye(length(obj.muNorm)),5);
            end
        end
        
        function value = fbeta(obj,x)
            if ~obj.wpNormal
                % Case beta distribution
                z = obj.prodtoforecast(x);
                %             value = betapdf(x/obj.factBeta,obj.a,obj.b);
                value = betapdf(z,obj.a,obj.b);
                value = prod(value);
            else
                % case function of a normal distribution
                xnorm = obj.hinv(x,'010');
                nx = length(xnorm);
                value = mvnpdf(xnorm,obj.muWPis,obj.sigWP);
            end
        end
        
        function value = fu(obj,x)
            if ~obj.flagU
                % case where we are not interested in the probability of u.
                % In this case, we set fu = 1.
                value = 1;
            else
                value = mvnpdf(x,obj.muU,obj.sigU);
            end
        end
        
        function value = fzeta(obj,Zbeta,Znorm)
            if obj.flagDecorr
                zeta = [Zbeta;Znorm];
                y = obj.hinv(zeta,'011');
                value = mvnpdf(y.');
            else
                value = obj.fbeta(Zbeta).*obj.fnorm(Znorm);
            end
        end
        
        function value = f(obj,lambda)
            u = lambda(obj.indU);
            Zbeta = lambda(obj.indBeta);
            Znorm = lambda(obj.indNorm);
            fofzeta = obj.fzeta(Zbeta,Znorm);
            value = fofzeta;
            if ~isempty(u)
                value = obj.fu(u)*fofzeta;
            end
        end
        
        function value = fload(obj,lambda)
            Zload = lambda(obj.indNorm);
            value = obj.fnorm(Zload);
        end
        
        %% Values of the derivatives for a given zeta
        function value = Dfnorm(obj,x)
            n = length(x);
            invSigNormX = obj.sigNorm\(x-obj.muNorm);
            %             value = - invSigNormX/sqrt((x-obj.muNorm).'*invSigNormX);
            value = 1/((2*pi)^(n/2)*sqrt(det(obj.sigNorm)))*(-invSigNormX)*exp(-1/2*(x-obj.muNorm).'*invSigNormX);
        end
        
        function value = Dfbeta(obj,x)
            if ~obj.wpNormal
                [z,dz_dx] = obj.prodtoforecast(x);
                % x = x/obj.factBeta;
                value = dz_dx.*1./beta(obj.a,obj.b).*((obj.a-1).*z.^(obj.a-2).*(1-z).^(obj.b-1)-(obj.b-1).*z.^(obj.a-1).*(1-z).^(obj.b-2));
            else
                % Case modelled by a function of a normal distribution
                xwp = obj.hinv(x,'010');
                n = length(xwp);
                invSigNormX = obj.sigWP\(xwp-obj.muWPis);
                % value = - invSigNormX/sqrt((x-obj.muNorm).'*invSigNormX);
                value = 1/((2*pi)^(n/2)*sqrt(det(obj.sigWP)))*(-invSigNormX)*exp(-1/2*(xwp-obj.muWPis).'*invSigNormX);
            end
        end
        
        function value = Dfu(obj,x)
            n = length(x);
            invSigNormX = obj.sigU\(x-obj.muU);
            value = 1/((2*pi)^(n/2)*sqrt(det(obj.sigU)))*(-invSigNormX)*exp(-1/2*(x-obj.muU).'*invSigNormX);
%             value = 1./(obj.sigU*sqrt(2*pi)).*(-(x-obj.muU)./obj.sigU).*exp(-(x-obj.muU).^2./(2*obj.sigU));
        end
        
        function value = Dfzeta(obj,Zbeta,Znorm)
            if obj.flagDecorr
                zeta = [Zbeta;Znorm];
                y = obj.hinv(zeta,'011');
                fy = mvnpdf(y.');
                value = -y*fy;
            else
                value = [obj.Dfbeta(Zbeta).*obj.fnorm(Znorm);...
                    obj.fbeta(Zbeta).*obj.Dfnorm(Znorm)];
            end
        end

        function value = Df(obj,lambda)
            u = lambda(obj.indU);
            Zbeta = lambda(obj.indBeta);
            Znorm = lambda(obj.indNorm);
            value = zeros(length(lambda),1);
            if ~isempty(u)
                value(obj.indU) = obj.Dfu(u)*obj.fbeta(Zbeta).*obj.fnorm(Znorm);
                value([obj.indBeta obj.indNorm]) = obj.fu(u)*obj.Dfzeta(Zbeta,Znorm);
%                 value(obj.indBeta) = obj.fu(u)*obj.Dfbeta(Zbeta).*obj.fnorm(Znorm);
%                 value(obj.indNorm) = obj.fu(u)*obj.fbeta(Zbeta).*obj.Dfnorm(Znorm);
            else
                value([obj.indBeta obj.indNorm]) = obj.Dfzeta(Zbeta,Znorm);
%                 value(obj.indBeta) = obj.Dfbeta(Zbeta).*obj.fnorm(Znorm);
%                 value(obj.indNorm) = obj.fbeta(Zbeta).*obj.Dfnorm(Znorm);
            end
        end
        
        %% Importance function
        function value = impFunc(obj,lambda)
            lambdat = obj.hinv(lambda);
            Zbeta = lambdat(obj.indBeta);
            Znorm = lambdat(obj.indNorm);
            x = [Zbeta;Znorm];
            invSigNormX = obj.sigX\(x-obj.muX);
            value = -(x-obj.muX).'*invSigNormX;
        end
        
        function value = diffImpFunc(obj,lambda)
            lambdat = obj.hinv(lambda);
            Zbeta = lambdat(obj.indBeta);
            Znorm = lambdat(obj.indNorm);
            x = [Zbeta;Znorm];
            invSigNormX = obj.sigX\(x-obj.muX);
            value = - 2*invSigNormX;
        end
        
        %% Going from random variables to beta distributed
        function lambda = h(obj,lambdatilde,varargin)
            % Lambda = h(lambdatilde)
            inc_u = 1;
            inc_wp = 1;
            inc_load = 1;
            if ~isempty(varargin)
                inLambda = varargin{1};
                inc_u = base2dec(inLambda(1),2);
                inc_load = base2dec(inLambda(3),2);
                inc_wp = base2dec(inLambda(2),2);
            end
            if inc_u == 0
                indWP = 1:obj.nwind;
                indLoad = obj.nwind+1:obj.nwind+obj.nload;
            else
                indWP = obj.indBeta;
                indLoad = obj.indNorm;
            end
            if obj.flagDecorr && inc_load && inc_wp
                % Then lambdatilde is actually in the y space, so we need
                % to transform it to the lambdatilde space first
                y = lambdatilde([indWP indLoad]);
                lambdatilde([indWP indLoad]) = obj.YtoX(y);
            end
            if obj.wpNormal
                % lambdatilde = lambda with normal random variables N(0,sigWP)
                % for WP
                % lambda = beta distributed ones
                lambda = lambdatilde;
                if ~isempty(indWP)
                    xwp = lambdatilde(indWP)-obj.muWPis;
                    u = normcdf(xwp,0,1);
                    u(u<=1e-9) = 0; % very small u:s are set to 0. See the hinv function.
                    if size(u,2) > 1
                        u = u.';
                    end
                    w = obj.WPcap.*betainv(u,obj.a,obj.b);
                    lambda(indWP) = w; % lambda = lambdatilde except for the new wp variables
                end
            else
                lambda = lambdatilde;
            end
        end
        
        function lambdatilde = hinv(obj,lambda,varargin)
            if obj.wpNormal
                % lambdatilde = lambda with normal random variables N(0,sigWP)
                % for WP
                % lambda = beta distributed ones
                inc_u = 1;
                inc_load = 1;
                inc_wp = 1;
                if ~isempty(varargin)
                    inLambda = varargin{1};
                    inc_u = base2dec(inLambda(1),2);
                    inc_load = base2dec(inLambda(3),2);
                    inc_wp = base2dec(inLambda(2),2);
                end
                
                if ~inc_u
                    % lambdatilde does not contain u
                    indWP = 1:obj.nwind;
                    indLoad = obj.nwind+1:obj.nwind+obj.nload;
                else
                    % lambdatilde contains also u
                    indWP = obj.indBeta;
                    indLoad = obj.indNorm;
                end
                lambdatilde = lambda;
                w = lambda(indWP);
                if ~isempty(w)
                    wnorm = w./obj.WPcap;
                    u = betacdf(wnorm,obj.a,obj.b);
                    % We do this to avoid numerical difficulties with the
                    % inverse function of the normal distribution when wind
                    % power does not produce anything.
                    u(u==0) = 1e-10;
                    xwp = norminv(u,0,1)+obj.muWPis;
                    lambdatilde(indWP) = xwp;
                else
                    xwp = [];
                end  
                if inc_wp == 1 && (inc_u+inc_load == 0) && ~isempty(xwp)
                    % only include wind power
                    lambdatilde = xwp;
                end
                if obj.flagDecorr
                    if inc_wp && inc_load
                        % Transform to the space of decorrelated variable Y
                        x = lambdatilde([indWP indLoad]);
                    elseif inc_wp
                        x = lambdatilde(indWP);
                    elseif inc_load
                        x = lambdatilde(indLoad);
                    end
                    lambdatilde = obj.XtoY(x,[num2str(inc_wp) num2str(inc_load)]);
                end
            else
                lambdatilde = lambda;
            end
        end
        
        function value = dh(obj,lambdatilde,varargin)
            if obj.wpNormal
                n = length(lambdatilde);
                s = ones(n,1);
                % for u or loads, lambda = lambdatilde, and the derivative is
                % one. We need only to change for the wind power
                lambda = obj.h(lambdatilde);
                if obj.flagDecorr
                    xwp0 = obj.XtoY(lambdatilde(obj.indBeta),'10');
                else
                    xwp0 = lambdatilde(obj.indBeta);
                end
                if ~isempty(varargin)
                    onlyWP = 1;
                else
                    onlyWP = 0;
                end
                if ~isempty(xwp0)
                    xwp = xwp0-obj.muWPis;
                    w = lambda(obj.indBeta);
                    wnorm = w./obj.WPcap;
                    fwnorm = betapdf(wnorm,obj.a,obj.b);
                    indZero = fwnorm == 0;
                    s(obj.indBeta) = obj.WPcap.*1./fwnorm.*normpdf(xwp,0,1);
                    s(indZero) = 0;
                end
                indi = 1:n;
                indj = 1:n;
                valuedh = sparse(indi,indj,s);
                if obj.flagDecorr
                    indZeta = [obj.indBeta obj.indNorm];
                    valuedh(indZeta,indZeta) = valuedh(indZeta,indZeta)*obj.Mchol;
                end
                if onlyWP
                    keyboard;
                    value = valuedh(obj.indBeta,obj.indBeta);
%                     value = sparse(indi,indj,s(obj.indBeta));
                else
                    value = valuedh;
                end
            else
                value = 1;
            end
        end
        
        function value = d2h(obj,lambdatilde,varargin)
            if obj.wpNormal
                n = length(lambdatilde);
                indwp = obj.indBeta;
                if size(indwp,2) > 1
                    indwp = indwp.';
                end
                % the only nonzero values are the ones for the wind power in
                % lambda
                lambda = obj.h(lambdatilde);
                if obj.flagDecorr
                    xwp0 = obj.XtoY(lambdatilde(obj.indBeta),'10');
                else
                    xwp0 = lambdatilde(obj.indBeta);
                end
                if ~isempty(varargin)
                    onlyWP = 1;
                else
                    onlyWP = 0;
                end
                if ~isempty(xwp0)
                    xwp = xwp0-obj.muWPis;
                    w = lambda(obj.indBeta);
                    wnorm = w./obj.WPcap;
                    fwnorm = betapdf(wnorm,obj.a,obj.b);
                    indZero = fwnorm == 0;
                    vals = obj.WPcap./fwnorm.*(-xwp.*normpdf(xwp,0,1)-obj.derivbetapdf(wnorm)./fwnorm.^2.*normpdf(xwp,0,1).^2);
                    vals(indZero) = 0;
                    subs = [indwp indwp indwp];
                    valued2h = sptensor(subs,vals,[n n n]);
                else
                    valued2h = 0;
                end
                if obj.flagDecorr
                    indZeta = [obj.indBeta obj.indNorm];
                    valued2hzeta = valued2h(indZeta,indZeta,indZeta);
                    valued2hzeta = ttm(valued2hzeta,{obj.Mchol.',obj.Mchol.'},[2 3]);
                    valued2h(indZeta,indZeta,indZeta) = valued2hzeta;
                end
                if onlyWP
                    keyboard;
                    nbeta = length(indwp);
                    subs = [(1:nbeta).' (1:nbeta).' (1:nbeta).'];
                    value = sptensor(subs,vals,[nbeta nbeta nbeta]);
                else
                    value = valued2h;
%                     subs = [indwp indwp indwp];
%                     value = sptensor(subs,vals,[n n n]);
                end
                
            else
                value = 0;
            end
        end
        
        function value = d3h(obj,lambdatilde,varargin)
            if obj.wpNormal
                n = length(lambdatilde);
                indwp = obj.indBeta;
                if size(indwp,2) > 1
                    indwp = indwp.';
                end
                % the only nonzero values are the ones for the wind power in
                % lambda
                
                lambda = obj.h(lambdatilde);
                if obj.flagDecorr
                    xwp0 = obj.XtoY(lambdatilde(obj.indBeta),'10');
                else
                    xwp0 = lambdatilde(obj.indBeta);
                end
                xwp = xwp0-obj.muWPis;
                w = lambda(obj.indBeta);
                wnorm = w./obj.WPcap;
                fwnorm = betapdf(wnorm,obj.a,obj.b);
                indZero = fwnorm == 0;
                fpwnorm = obj.derivbetapdf(wnorm);
                phi = normpdf(xwp,0,1);
                phip = -xwp.*normpdf(xwp,0,1);
                fppwnorm = obj.deriv2betapdf(wnorm);
                phipp = (xwp.^2-1).*normpdf(xwp,0,1);
                term1 = -fpwnorm.*phi./(fwnorm.^3).*(phip-fpwnorm.*phi.^2./(fwnorm.^2));
                term2 = 1./fwnorm.*(phipp-((fppwnorm.*fwnorm-2*fpwnorm.^2)./(fwnorm.^4).*phi.^3+...
                    2*fpwnorm.*phip.*phi./(fwnorm.^2)));
                vals = obj.WPcap.*(term1+term2);
                vals(indZero) = 0;
                if ~isempty(varargin)
                    onlyWP = 1;
                else
                    onlyWP = 0;
                end
                subs = [indwp indwp indwp indwp];
                valued3h = sptensor(subs,vals,[n n n n]);
                if obj.flagDecorr
                    indZeta = [obj.indBeta obj.indNorm];
                    valued3hzeta = valued3h(indZeta,indZeta,indZeta,indZeta);
                    valued3hzeta = ttm(valued3hzeta,{obj.Mchol.',obj.Mchol.',obj.Mchol.'},[2 3 4]);
                    valued3h(indZeta,indZeta,indZeta,indZeta) = valued3hzeta;
                end
                if onlyWP
                    keyboard;
                    nbeta = length(indwp);
                    subs = [(1:nbeta).' (1:nbeta).' (1:nbeta).' (1:nbeta).'];
                    value = sptensor(subs,vals,[nbeta nbeta nbeta nbeta]);
                else
                    value = valued3h;
                end
                
            else
                value = 0;
            end
        end
        
        function value = derivbetapdf(obj,x)
            % Derivative of the beta pdf with the parameters in the object
            value = 1./beta(obj.a,obj.b).*((obj.a-1).*x.^(obj.a-2).*(1-x).^(obj.b-1)-(obj.b-1).*x.^(obj.a-1).*(1-x).^(obj.b-2));
            indOut = x<0 | x>1;
            if nnz(indOut) ~= 0
                % there shouldnt be any x outside [0,1]. Firewall here to
                % verify that
                keyboard;
            end
            value(indOut) = 0;
        end
        
        function value = deriv2betapdf(obj,x)
            % Derivative of the beta pdf with the parameters in the object
            value = 1./beta(obj.a,obj.b).*((obj.a-1).*(obj.a-2).*x.^(obj.a-3).*(1-x).^(obj.b-1) ...
                -2*(obj.b-1).*(obj.a-1).*x.^(obj.a-2).*(1-x).^(obj.b-2)...
                + (obj.b-1).*(obj.b-2).*x.^(obj.a-1).*(1-x).^(obj.b-3));
            indOut = x<0 | x>1;
            if nnz(indOut) ~= 0
                % there shouldnt be any x outside [0,1]. Firewall here to
                % verify that
                keyboard;
            end
            value(indOut) = 0;
        end
        
        %% Function for transforming from tilde(zeta) = X to Y by
        % X = mu + M*Y, Y ~ N(0,I) and M M^T = Sigma.
        
        function values = YtoX(obj,valuesY,varargin)
            % X = mu + M*Y
            inc_load = 1;
            inc_wp = 1;
            if ~isempty(varargin)
                inLambda = varargin{1};
                inc_wp = base2dec(inLambda(1),2);
                inc_load = base2dec(inLambda(2),2);
            end
            if inc_wp && inc_load
                values = bsxfun(@plus,obj.muX,obj.Mchol*valuesY);
            elseif inc_wp
                values = bsxfun(@plus,obj.muXwp,obj.Mwp*valuesY);
            elseif inc_load
                values = bsxfun(@plus,obj.muNorm,obj.Ml*valuesY);
            end
        end
        
        function values = XtoY(obj,valuesX,varargin)
            % Y = M^(-1) (X-mu)
            inc_load = 1;
            inc_wp = 1;
            if ~isempty(varargin)
                inLambda = varargin{1};
                inc_wp = base2dec(inLambda(1),2);
                inc_load = base2dec(inLambda(2),2);
            end
            
            if inc_wp && inc_load
                Minv = inv(obj.Mchol);
                values = Minv*bsxfun(@minus,valuesX,obj.muX);
            elseif inc_wp 
                Minv = inv(obj.Mwp);
                values = Minv*bsxfun(@minus,valuesX,obj.muXwp);
            elseif inc_load
                Minv = inv(obj.Ml);
                values = Minv*bsxfun(@minus,valuesX,obj.muNorm);
            end
                
        end
        
        function values = DyX(obj)
            values = obj.Mchol;
        end
        
        %% Generate forecast
        function samples = getSamples(obj,n)
            if ~obj.modeStudentt
                samplesNorm = mvnrnd(obj.muNorm,obj.sigNorm,n);
            else
                nload = length(obj.muNorm);
                samplestudent = mvtrnd(obj.sigNorm,5,n);
                samplesNorm = zeros(n,nload);
                for i = 1:n
                    samplesNorm(i,:) = (samplestudent(i,:).'+obj.muNorm);
                end
            end
            if ~obj.wpNormal
                nbWind = length(obj.a);
                a_array = repmat(obj.a.',n,1);
                b_array = repmat(obj.b.',n,1);
                samplesBeta = betarnd(a_array,b_array,n,nbWind);
                for i = 1:n
                    samplesBeta(i,:) = (obj.WPcap.').*samplesBeta(i,:);
                end
            else
                nbeta = length(obj.indBeta);
                samplesBeta = mvnrnd(obj.muWPis,obj.sigWP,n);
            end
            samples = [samplesBeta samplesNorm].';
        end
        
        %% Draw distributions
        function [hfig,xvalues,probs] = draw(obj,option)
            % option = 1 draws the distribution for the beta variables
            %   and the associated normally distributed variables if they
            %   exist
            % option = 2 draws the distribution for the loads
            nvar = length(obj.indBeta);
            nbPts = 100;
            if option == 2
                xload = linspace(obj.muNorm-3*sqrt(obj.sigNorm),obj.muNorm+3*sqrt(obj.sigNorm),nbPts).';
                xvalues = xload;
                probs = normpdf(xload,obj.muNorm,sqrt(obj.sigNorm));
                hfig = figure;
                plot(xload,probs,'k');
                xlabel('P_{load}');
                ylabel('Density');
            elseif option == 1
                % First draw the beta distribution
                xbeta = linspace(0,1,nbPts);
                probbeta = zeros(nbPts,1);
                for i = 1:nbPts
                    probbeta(i) = betapdf(xbeta(i),obj.a,obj.b);
                end
                wp_beta = xbeta*obj.WPcap;
                xvalues = wp_beta;
                
                if obj.wpNormal
                    xnorm = linspace(-5,5,nbPts);
                    xnorm = obj.hinv(wp_beta,'010');
                    xwp = zeros(nvar,nbPts);
                    probs = zeros(nbPts,1);
                    for j = 1:nbPts
                        xwp(:,j) = obj.h(xnorm(j)*ones(nvar,1));
                        mu = zeros(nvar,1);
                        probs(j) = mvnpdf(xnorm(j)*ones(nvar,1),mu,obj.sigWP);
                    end
                end
                
                hfig = figure;
                ax1 = gca;
                hold on
                plot(wp_beta,probbeta,'k');
                xlabel('Wind power production P_{g2}');
                ylabel('Density');
                xlim([0 obj.WPcap]);
                ax2 = axes('Position',get(ax1,'Position'),...
                    'XAxisLocation','top',...
                    'Color','none',...
                    'XColor','k','YColor','k');
                xlabel(ax2,'x_{g2}');
                % Second axis is from 0 to 1
                % First axis from 0 to WPcap
                tickVal_tilde = [-10;-5;0;5];
                tickVal_norm = zeros(length(tickVal_tilde),1);
                for i = 1:length(tickVal_tilde)
                    tickVal_norm(i) = obj.h(tickVal_tilde(i),1)/obj.WPcap; % corresponding values between 0 and 1
                end
                set(ax2,'XTick',tickVal_norm);
                set(ax2,'XTickLabel',tickVal_tilde);
                linkaxes([ax1 ax2],'y');
            end
        end
        
        function writeStochSettings(obj,scenarName)
            % This function write the current stochastic parameters in the
            % file scenarName
            fid = fopen(sprintg('%s.txt',scenarName),'w');
            fprintf(fid,'\nStochastic model\n');
            fprintf(fid,'--------\n');
            fprintf(fid,'\nModel for wind power:  ');
            if obj.wpNormal
                fprintf(fid,'Gaussian (parameters zeta tilde)');
            else
                fprintf(fid,'non Gaussian (parameters zeta)');
            end
            
            fprintf(fid,'\n*Wind power: Alpha, beta, sigWP*\n');
            fprintf(fid,'-------------------------------------------------------------------------\n');
            fprintf(fid,'Alpha  :  ');fprintf(fid,'%6d\t',obj.a);fprintf(fid,'\n');
            fprintf(fid,'Beta   :  ');fprintf(fid,'%6d\t',obj.b);fprintf(fid,'\n');
            fprintf(fid,'Covariance matrix:  ');
            for i=1:size(obj.sigWP,1)
                for j=1:size(obj.sigWP,2)
                    fprintf(fid,'%6d\t',obj.sigWP(i,j));
                end
                fprintf(fid,'\n');
            end
            fprintf(fid,'-------------------------------------------------------------------------\n');
            
            fprintf(fid,'*Average load and covariance matrix*\n');
            fprintf(fid,'-------------------------------------------------------------------------\n');
            fprintf(fid,'Averages:  ');fprintf(fid,'%6d\t',obj.muNorm);fprintf(fid,'\n');
            fprintf(fid,'Covariance matrix:  ');
            for i=1:size(obj.sigNorm,1)
                for j=1:size(obj.sigNorm,2)
                    fprintf(fid,'%6d\t',obj.sigNorm(i,j));
                end
                fprintf(fid,'\n');
            end
            fprintf(fid,'\n');
        end
        
        function hfig = drawXvsPg(obj)
            nbPts = 50;
            Pg = linspace(0,obj.WPcap,nbPts);
            xs = zeros(nbPts,1);
            for i = 1:nbPts
                xs(i) = obj.hinv(Pg(i),'010');
            end
            hfig = figure;
            plot(Pg,xs,'k');
        end
    end
end