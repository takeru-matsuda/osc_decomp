function [osc_param,osc_AIC,osc_mean,osc_cov,osc_phase] = osc_decomp_uni(y,fs,MAX_OSC,MAX_AR,algorithm,grad)
%
% Input:
%    y:             univariate time series (1 times T)
%    fs:            sampling frequency (scalar)
%    MAX_OSC:       maximum number of oscillation components (scalar)
%    MAX_AR:        maximum AR degree (scalar), default = 2*MAX_OSC
%    algorithm:     optimization algorithm (string), default = 'quasi-Newton'
%    grad:          gradient computation method in optimization (boolean), default = false
%                   true -> Kalman filter, false -> numerical gradient
%
% Output:
%    osc_param:     estimated parameters of the oscillator model (MAX_OSC times 3*MAX_OSC+1)
%                   osc_param(K,1:K) is the estimated a_1,...,a_K of the oscillator model with K oscillation components
%                   osc_param(K,K+1:2*K) is the estimated f_1,...f_K of the oscillator model with K oscillation components
%                   osc_param(K,2*K+1:3*K) is the estimated sigma_1^2,...,sigma_K^2 of the oscillator model with K oscillation components
%                   osc_param(K,3*K+1) is the estimated tau^2 of the oscillator model with K oscillation components
%    osc_AIC:       AIC of the oscillator model (1 times MAX_OSC)
%                   AIC_osc(K) is the AIC of the oscillator model with K oscillation components
%    osc_mean:      smoothed coordinate of each oscillator (2*MAX_OSC times T times MAX_OSC)
%                   decomp_mu(2*k-1:2*k,:,K) is the smoothed coordinate of the k-th oscillator in the decomposition into K components
%    osc_cov:       smoothed covariance of each oscillator (2*MAX_OSC times 2*MAX_OSC times T times MAX_OSC)
%                   decomp_cov(2*k-1:2*k,2*k-1:2*k,:,K) is the smoothed covariance of the k-th oscillator in the decomposition into K components
%    osc_phase:     estimated phase of each oscillation component (MAX_OSC times T times MAX_OSC)
%                   phi_prop(k,:,K) is the estimated phase of the k-th oscillator in the decomposition into K components
%

    arguments
        y
        fs (1,1)
        MAX_OSC (1,1) {mustBeInteger(MAX_OSC),mustBePositive(MAX_OSC)} = 5
        MAX_AR (1,1) {mustBeInteger(MAX_AR),mustBePositive(MAX_AR)} = max(20,2*MAX_OSC)
        algorithm = 'quasi-newton'
        grad = false
    end
    T = length(y);
    if nargin == 3
        MAX_AR = 2*MAX_OSC;
    else
        MAX_AR = max(MAX_AR,2*MAX_OSC);
    end
    [ARwithnoise_param,ARwithnoise_AIC] = AR_fit(y,MAX_AR);
    osc_param = zeros(MAX_OSC,3*MAX_OSC+1);
    osc_AIC = zeros(1,MAX_OSC);
    osc_mean = zeros(2*MAX_OSC,T,MAX_OSC);
    osc_cov = zeros(2*MAX_OSC,2*MAX_OSC,T,MAX_OSC);
    osc_phase = zeros(MAX_OSC,T,MAX_OSC);
	for K=1:MAX_OSC
        K
        [~,ARdeg] = min(ARwithnoise_AIC);
        tmp = roots([1 ARwithnoise_param(ARdeg,1:ARdeg)]);
        if nnz(imag(tmp)>=0) >= K
            z0 = tmp;
            E0 = ARwithnoise_param(ARdeg,ARdeg+1);
            R0 = ARwithnoise_param(ARdeg,ARdeg+2);
            optARdeg = ARdeg;
        else
            minAIC = inf;
            minAIC2 = inf;
            for ARdeg=K:2*K
                tmp = roots([1 ARwithnoise_param(ARdeg,1:ARdeg)]);
                if ARdeg-nnz(imag(tmp))/2 == K && ARwithnoise_AIC(ARdeg) < minAIC
                    z0 = tmp;
                    E0 = ARwithnoise_param(ARdeg,ARdeg+1);
                    R0 = ARwithnoise_param(ARdeg,ARdeg+2);
                    minAIC = ARwithnoise_AIC(ARdeg);
                    optARdeg = ARdeg;
                end
                if ARdeg-nnz(imag(tmp))/2 >= K && ARwithnoise_AIC(ARdeg) < minAIC2
                    z1 = tmp;
                    E1 = ARwithnoise_param(ARdeg,ARdeg+1);
                    R1 = ARwithnoise_param(ARdeg,ARdeg+2);
                    minAIC2 = ARwithnoise_AIC(ARdeg);
                    optARdeg2 = ARdeg;
                end
            end
            if minAIC == inf
                if minAIC2 == inf
                    warning('no AR model with %d oscillators',K);
                end
                z0 = z1;
                E0 = E1;
                R0 = R1;
                optARdeg = optARdeg2;
                [~,I] = sort(abs(z0),'descend');
                z0 = z0(I);
            end
        end
        VV = zeros(optARdeg,optARdeg);
        for j=1:optARdeg
            for i=1:optARdeg
                VV(i,j) = z0(j)^(1-i);
            end
        end
        QQ = inv(VV)*[E0 zeros(1,optARdeg-1); zeros(optARdeg-1,optARdeg)]*inv(VV)';
        [~,I] = sort(diag(real(QQ))./(1-abs(z0).^2),'descend');
        z0 = z0(I);
        init_a = zeros(1,K);
        init_theta = zeros(1,K);
        kk = 1;
        for k=1:K
            init_a(k) = abs(z0(kk));
            init_theta(k) = abs(angle(z0(kk)));
            if imag(z0(kk)) == 0
                kk = kk+1;
            else
                kk = kk+2;
            end
        end
        if mod(T,2) == 0
        	freq = [2*pi/T*(0:T/2-1) pi];
	    else
    	    freq = [2*pi/T*(0:(T-1)/2)];
        end
        P = zeros(length(freq),K+1);
        for k=1:K
            a = init_a(k);
            theta = init_theta(k);
            A = (1-2*a^2*cos(theta)^2+a^4*cos(2*theta))/a/(a^2-1)/cos(theta);
            b = (A-2*a*cos(theta)+sign(cos(theta))*sqrt((A-2*a*cos(theta))^2-4))/2;
            for j=1:length(freq)
                P(j,k) = -a*cos(theta)/b*abs(1+b*exp(-1i*freq(j)))^2/abs(1-2*a*cos(theta)*exp(-1i*freq(j))+a^2*exp(-2*1i*freq(j))).^2/2/pi;
            end
        end
        P(:,K+1) = 1/2/pi;
        p = zeros(length(freq),1);
        for j=1:length(freq)
            p(j) = abs(y*exp(-1i*freq(j)*(0:T-1)'))^2/2/pi/T;
        end
        weight = whittle_uni_fit(P,p);
        init_sigma2 = weight(1:K)';
        init_tau2 = weight(K+1);

        if nargin >= 5
            switch lower(algorithm)
                case 'quasi-newton'
                    if nargin == 6 && grad == true
                        options = optimoptions('fminunc','Algorithm','quasi-newton','MaxFunEvals',10000,'SpecifyObjectiveGradient',true);
                        param = fminunc(@(param)osc_uni_prof_ll(y,param,init_theta,true),[atanh(2*init_a-1) zeros(1,K) log(init_sigma2/init_tau2)],options);
                    else
                        options = optimoptions('fminunc','Algorithm','quasi-newton','MaxFunEvals',10000,'SpecifyObjectiveGradient',false);
                        param = fminunc(@(param)osc_uni_prof_ll(y,param,init_theta,false),[atanh(2*init_a-1) zeros(1,K) log(init_sigma2/init_tau2)],options);
                    end
                case 'trust-region'
                    options = optimoptions('fminunc','Algorithm','trust-region','MaxFunEvals',10000,'SpecifyObjectiveGradient',true);
                    param = fminunc(@(param)osc_uni_prof_ll(y,param,init_theta,true),[atanh(2*init_a-1) zeros(1,K) log(init_sigma2/init_tau2)],options);
                case 'cg'
                    param = minimize([atanh(2*init_a-1) zeros(1,K) log(init_sigma2/init_tau2)]',@(param)osc_uni_prof_ll(y,param,init_theta,true),-10000)';
                case 'em'
                    param = osc_em(y,[init_a init_theta init_sigma2 init_tau2]);
                otherwise
                    error('osc_decomp_multi: unknown algorithm')
            end
        else
            options = optimoptions('fminunc','Algorithm','quasi-newton','MaxFunEvals',10000,'SpecifyObjectiveGradient',false);
            param = fminunc(@(param)osc_uni_prof_ll(y,param,init_theta,false),[atanh(2*init_a-1) zeros(1,K) log(init_sigma2/init_tau2)],options);
        end
    	[mll,~,osc_tau2] = osc_uni_prof_ll(y,param,init_theta,false);
    	osc_AIC(K) = 2*mll+2*(3*K+1);

        param(K+1:2*K) = init_theta+tanh(param(K+1:2*K))*pi;
	    [~,I] = sort(abs(param(K+1:2*K)));
    	osc_a = (tanh(param(I))+1)/2;
    	osc_f = abs(param(K+I))*fs/2/pi;
    	osc_sigma2 = exp(param(2*K+I))*osc_tau2;
		osc_param(K,1:3*K+1) = [osc_a osc_f osc_sigma2 osc_tau2];
        [osc_mean(1:2*K,:,K),osc_cov(1:2*K,1:2*K,:,K)] = osc_smooth(y,fs,osc_param(K,1:3*K+1));
	    for k=1:K
        	osc_phase(k,:,K) = atan2(osc_mean(2*k,:,K),osc_mean(2*k-1,:,K));
    	end
	end
end


