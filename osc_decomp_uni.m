function [osc_param,osc_AIC,osc_mean,osc_cov,osc_phase] = osc_decomp_uni(y,fs,MAX_OSC,algorithm,obj_func,ss_interval)
% Input:
%    y:             Univariate time series (1 x T)
%    fs:            Sampling frequency (scalar)
%    MAX_OSC:       Maximum number of oscillators (integer), default = 5
%    algorithm:     Optimization algorithm (string), default = 'quasi-newton_grad0'
%                   Supported: 'quasi-newton_grad0', 'quasi-newton_grad1', 'trust-region', 'cg', 'adam', 'em'
%    obj_func:      Objective function (string), default = 'whittle'
%                   Supported: 'kalman', 'whittle'
%    ss_interval:   Subsampling interval for frequency domain (integer), default = 1
%                   (Only valid when obj_func = 'whittle')
%
% Output:
%    osc_param:     Estimated parameters (MAX_OSC x 3*MAX_OSC+1)
%                   [a_1,...,a_K, f_1,...,f_K, sigma2_1,...,sigma2_K, tau2]
%    osc_AIC:       AIC of the oscillator model for each K (1 x MAX_OSC)
%    osc_mean:      Smoothed state means (2*MAX_OSC x T x MAX_OSC)
%    osc_cov:       Smoothed state covariances (2*MAX_OSC x 2*MAX_OSC x T x MAX_OSC)
%    osc_phase:     Estimated instantaneous phases (MAX_OSC x T x MAX_OSC)

    arguments
        y
        fs (1,1)
        MAX_OSC (1,1) {mustBeInteger(MAX_OSC),mustBePositive(MAX_OSC)} = 5
        algorithm = 'quasi-newton_grad1'
        obj_func = 'whittle'
        ss_interval = 1
    end
    
    T = length(y);
    if T<500
        obj_func = 'kalman';
    end

    MAX_AR = max(10,2*MAX_OSC);
    [ARwithnoise_param,ARwithnoise_AIC] = AR_fit(y,MAX_AR);
    
    osc_param = zeros(MAX_OSC,3*MAX_OSC+1);
    osc_AIC = zeros(1,MAX_OSC);
    osc_mean = zeros(2*MAX_OSC,T,MAX_OSC);
    osc_cov = zeros(2*MAX_OSC,2*MAX_OSC,T,MAX_OSC);
    osc_phase = zeros(MAX_OSC,T,MAX_OSC);
    
	for K = 1:MAX_OSC
        K
        [~,ARdeg] = min(ARwithnoise_AIC);
        tmp = roots([1 ARwithnoise_param(ARdeg,1:ARdeg)]);
        if nnz(imag(tmp)>=0) >= K
            z0 = tmp; E0 = ARwithnoise_param(ARdeg,ARdeg+1); R0 = ARwithnoise_param(ARdeg,ARdeg+2); optARdeg = ARdeg;
        else
            minAIC = inf; minAIC2 = inf;
            for ARdeg = K:2*K
                tmp = roots([1 ARwithnoise_param(ARdeg,1:ARdeg)]);
                if ARdeg-nnz(imag(tmp))/2 == K && ARwithnoise_AIC(ARdeg) < minAIC
                    z0 = tmp; E0 = ARwithnoise_param(ARdeg,ARdeg+1); R0 = ARwithnoise_param(ARdeg,ARdeg+2);
                    minAIC = ARwithnoise_AIC(ARdeg); optARdeg = ARdeg;
                end
                if ARdeg-nnz(imag(tmp))/2 >= K && ARwithnoise_AIC(ARdeg) < minAIC2
                    z1 = tmp; E1 = ARwithnoise_param(ARdeg,ARdeg+1); R1 = ARwithnoise_param(ARdeg,ARdeg+2);
                    minAIC2 = ARwithnoise_AIC(ARdeg); optARdeg2 = ARdeg;
                end
            end
            if minAIC == inf
                if minAIC2 == inf, warning('No AR model found with %d oscillators.', K); end
                z0 = z1; E0 = E1; R0 = R1; optARdeg = optARdeg2;
                [~,I] = sort(abs(z0),'descend'); z0 = z0(I);
            end
        end
        
        VV = zeros(optARdeg,optARdeg);
        for j=1:optARdeg, for i=1:optARdeg, VV(i,j) = z0(j)^(1-i); end, end
        QQ = inv(VV)*[E0 zeros(1,optARdeg-1); zeros(optARdeg-1,optARdeg)]*inv(VV)';
        [~,I] = sort(diag(real(QQ))./(1-abs(z0).^2),'descend');
        z0 = z0(I);
        
        init_a = zeros(1,K); init_theta = zeros(1,K); kk = 1;
        for k=1:K
            init_a(k) = abs(z0(kk)); init_theta(k) = abs(angle(z0(kk)));
            if imag(z0(kk)) == 0, kk = kk+1; else, kk = kk+2; end
        end
        
        m = floor(T/2) + 1;
        freq0 = 2*pi*(0:m-1)'/T;
        P = zeros(length(freq0),K+1);
        for k=1:K
            a = init_a(k); theta = init_theta(k);
            A = (1-2*a^2*cos(theta)^2+a^4*cos(2*theta))/a/(a^2-1)/cos(theta);
            b = (A-2*a*cos(theta)+sign(cos(theta))*sqrt((A-2*a*cos(theta))^2-4))/2;
            for j=1:length(freq0)
                P(j,k) = -a*cos(theta)/b*abs(1+b*exp(-1i*freq0(j)))^2/abs(1-2*a*cos(theta)*exp(-1i*freq0(j))+a^2*exp(-2*1i*freq0(j))).^2/2/pi;
            end
        end
        P(:,K+1) = 1/2/pi;
        Y_fft = fft(y);
        p_full = abs(Y_fft(1:m)).^2/T;
        options = optimoptions('fminunc','Algorithm','quasi-newton','SpecifyObjectiveGradient',true);
        beta = fminunc(@(b)whittle_regression_uni(b, P, p_full(:)), zeros(K+1,1), options);
        weight = exp(beta);
        init_sigma2 = weight(1:K)';
        init_tau2 = weight(K+1);
        
        param0 = [atanh(2*init_a-1) zeros(1,K) log(init_sigma2/init_tau2)]';

        switch lower(obj_func)
            case 'kalman'
                obj = @(param)kalman_uni_ll_wrapper(param,init_theta,y);

            case 'whittle'
                freq_sub = freq0(1:ss_interval:end);
                p_sub = p_full(1:ss_interval:end);
                obj = @(param)whittle_uni_ll_wrapper(param, init_theta, freq_sub, p_sub);

            otherwise
                error('osc_decomp_uni: Unknown obj_func specified.')
        end
        
        switch lower(algorithm)
            case 'quasi-newton_grad0'
                options = optimoptions('fminunc','Algorithm','quasi-newton','MaxFunEvals',1000*K,'SpecifyObjectiveGradient',false,...
                                       'FunctionTolerance',1e-4,'StepTolerance',1e-4,'OptimalityTolerance',1e-3);
                param = fminunc(obj, param0, options);

            case 'quasi-newton_grad1'
                options = optimoptions('fminunc','Algorithm','quasi-newton','MaxFunEvals',1000*K,'SpecifyObjectiveGradient',true,...
                                       'FunctionTolerance',1e-4,'StepTolerance',1e-4,'OptimalityTolerance',1e-3);
                param = fminunc(obj, param0, options);

            case 'trust-region'
                options = optimoptions('fminunc','Algorithm','trust-region','MaxFunEvals',10000,'SpecifyObjectiveGradient',true);
                param = fminunc(obj, param0, options);
                
            case 'cg'
                param = minimize_noshow(param0,obj,-1000);
                
            case 'adam'
                param = param0;
                N_freq = length(freq0);
                subsample_size = max(10, ceil(N_freq / ss_interval));
                
                maxIters = 2000;
                minIters = 200;
                tol = 1e-5;
                
                lr = 0.01; 
                beta1 = 0.9; beta2 = 0.999; eps_adam = 1e-8;
                
                m_adam = zeros(size(param));
                v_adam = zeros(size(param));
                
                mll_history = zeros(maxIters, 1);
                
                for t = 1:maxIters
                    idx = randi([2, N_freq], subsample_size, 1);
                    freq_batch = freq0(idx);
                    p_batch = p_full(idx);
                    
                    [mll_t, ~, grad] = whittle_uni_ll(param, init_theta, freq_batch, p_batch);
                    mll_history(t) = mll_t;
                    
                    if norm(grad) > 5.0
                        grad = grad * (5.0 / norm(grad));
                    end
                    
                    m_adam = beta1 * m_adam + (1 - beta1) * grad;
                    v_adam = beta2 * v_adam + (1 - beta2) * (grad.^2);
                    
                    mhat = m_adam / (1 - beta1^t);
                    vhat = v_adam / (1 - beta2^t);
                    
                    param = param - lr * mhat ./ (sqrt(vhat) + eps_adam);
                    
                    if t >= minIters && mod(t, 10) == 0
                        past_mean = mean(mll_history(t-19:t-10));
                        curr_mean = mean(mll_history(t-9:t));
                        if abs(past_mean - curr_mean) < tol
                            break;
                        end
                    end
                end
                fprintf('Adam converged in %d iterations',t)

            case 'em'
                param = osc_em(y,[init_a init_theta init_sigma2 init_tau2]);
                
            otherwise
                error('osc_decomp_uni: Unknown algorithm specified.')
        end
        
        param = param(:)';
    	[mll,osc_tau2] = kalman_uni_ll(param,init_theta,y);
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

function [obj, grad] = whittle_regression_uni(b, X, y)
    exp_b = exp(b);
    sp = X * exp_b;
    obj = sum(log(sp) + y ./ sp);    
    if nargout > 1
        weight = (1 ./ sp) - (y ./ (sp.^2));
        grad = (X' * weight) .* exp_b;
    end
end

function [mll, grad] = kalman_uni_ll_wrapper(param, init_theta, y)
    if nargout > 1
        [mll, ~, grad] = kalman_uni_ll(param, init_theta, y);
    else
        mll = kalman_uni_ll(param, init_theta, y);
    end
end

function [mll, grad] = whittle_uni_ll_wrapper(param, init_theta, freq, p)
    if nargout > 1
        [mll, ~, grad] = whittle_uni_ll(param, init_theta, freq, p);
    else
        mll = whittle_uni_ll(param, init_theta, freq, p);
    end
end
