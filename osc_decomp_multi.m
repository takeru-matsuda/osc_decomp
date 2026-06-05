function [osc_param,osc_AIC,osc_mean,osc_cov,osc_phase] = osc_decomp_multi(Y,fs,MAX_OSC,algorithm,obj_func,ss_interval)
% Input:
%    Y:             Multivariate time series (J x T)
%    fs:            Sampling frequency (scalar)
%    MAX_OSC:       Maximum number of oscillators (integer), default = 5
%    algorithm:     Optimization algorithm (string), default = 'quasi-newton_grad0'
%                   Supported: 'quasi-newton_grad0', 'quasi-newton_grad1', 'trust-region', 'adam'
%    obj_func:      Objective function (string), default = 'whittle'
%                   Supported: 'kalman', 'whittle'
%    ss_interval:   Subsampling interval for frequency domain (integer), default = 1
%                   (Only valid when obj_func = 'whittle')
%
% Output:
%    osc_param:     estimated parameters of the oscillator model (MAX_OSC x (2*J+1)*MAX_OSC+1)
%                   osc_param(K,1:K) is the estimated a_1,...,a_K of the oscillator model with K oscillation components
%                   osc_param(K,K+1:2*K) is the estimated f_1,...f_K of the oscillator model with K oscillation components
%                   osc_param(K,2*K+1:3*K) is the estimated sigma_1^2,...,sigma_K^2 of the oscillator model with K oscillation components
%                   osc_param(K,3*K+1:(2*J+1)*K) is the estimated c_{21,1},c_{21,2},...,c_{J1,1},c_{J1,2},...,c_{JK,1},c_{JK,2} of the oscillator model with K oscillation components
%                   osc_param(K,(2*J+1)*K+1) is the estimated tau^2 of the oscillator model with K oscillation components
%    osc_AIC:       AIC of the oscillator model (1 x MAX_OSC)
%                   AIC_osc(K) is the AIC of the oscillator model with K oscillation components
%    osc_mean:      smoothed coordinate of each oscillator (2*MAX_OSC x T x MAX_OSC)
%                   decomp_mu(2*k-1:2*k,:,K) is the smoothed coordinate of the k-th oscillator in the decomposition into K components
%    osc_cov:       smoothed covariance of each oscillator (2*MAX_OSC x 2*MAX_OSC x T x MAX_OSC)
%                   decomp_cov(2*k-1:2*k,2*k-1:2*k,:,K) is the smoothed covariance of the k-th oscillator in the decomposition into K components
%    osc_phase:     estimated phase of each oscillation component (MAX_OSC x T x MAX_OSC)
%                   phi_prop(k,:,K) is the estimated phase of the k-th oscillator in the decomposition into K components
%
    arguments
        Y
        fs (1,1)
        MAX_OSC (1,1) {mustBeInteger(MAX_OSC),mustBePositive(MAX_OSC)} = 5
        algorithm = 'quasi-newton_grad1'
        obj_func = 'whittle'
        ss_interval = 1
    end
    
	J = size(Y, 1);
	T = size(Y, 2);
    if T<200
        obj_func = 'kalman';
    end
    
    MAX_VAR = max(10, ceil(2 * MAX_OSC / J));
    [VARwithnoise_A, VARwithnoise_E, VARwithnoise_r, VARwithnoise_AIC] = VAR_fit(Y, MAX_VAR);
    
    osc_param = zeros(MAX_OSC, (2*J+1)*MAX_OSC+1);
    osc_AIC   = zeros(1, MAX_OSC);
    osc_mean  = zeros(2*MAX_OSC, T, MAX_OSC);
    osc_cov   = zeros(2*MAX_OSC, 2*MAX_OSC, T, MAX_OSC);
    osc_phase = zeros(MAX_OSC, T, MAX_OSC);
    
	for K = 1:MAX_OSC
        K
    	[~, ARdeg] = min(VARwithnoise_AIC);
	    [Vtmp, tmp] = polyeig_VAR(VARwithnoise_A(:, 1:J*ARdeg, ARdeg));
        
    	if nnz(imag(tmp) >= 0) >= K
	        V0 = Vtmp; z0 = tmp;
        	E0 = VARwithnoise_E(:,:,ARdeg); R0 = VARwithnoise_r(ARdeg);
        	minAIC = VARwithnoise_AIC(ARdeg); optARdeg = ARdeg;
        else
    	    minAIC = inf; minAIC2 = inf; minK = inf;
        	for ARdeg = ceil(K/J):MAX_VAR
                [Vtmp, tmp] = polyeig_VAR(VARwithnoise_A(:, 1:J*ARdeg, ARdeg));
                n_osc = J*ARdeg - nnz(imag(tmp))/2;
                
        	    if n_osc == K && VARwithnoise_AIC(ARdeg) < minAIC
            	    V0 = Vtmp; z0 = tmp;
                	E0 = VARwithnoise_E(:,:,ARdeg); R0 = VARwithnoise_r(ARdeg);
                	minAIC = VARwithnoise_AIC(ARdeg); optARdeg = ARdeg;
                end
    	        if n_osc > K && (n_osc < minK || (n_osc == minK && VARwithnoise_AIC(ARdeg) < minAIC2))
                    V1 = Vtmp; z1 = tmp;
                    E1 = VARwithnoise_E(:,:,ARdeg); R1 = VARwithnoise_r(ARdeg);
                	minAIC2 = VARwithnoise_AIC(ARdeg); optARdeg2 = ARdeg;
        	        minK = n_osc;
                end
            end
        	if minAIC == inf
            	if minAIC2 == inf, warning('No VAR model found with %d oscillators.', K); end
        	    V0 = V1; z0 = z1; E0 = E1; R0 = R1; optARdeg = optARdeg2;
            end
        end
        
	    VV = zeros(J*optARdeg, J*optARdeg);
    	for j = 1:J*optARdeg
        	for i = 1:optARdeg
            	VV((i-1)*J+1:i*J, j) = z0(j)^(1-i) * V0(:,j);
    	    end
	    end
    	QQ = inv(VV) * [E0 zeros(J, J*(optARdeg-1)); zeros(J*(optARdeg-1), J*optARdeg)] * inv(VV)';
	    [~, I] = sort(diag(real(QQ)) ./ (1 - abs(z0).^2), 'descend');
    	V0 = V0(:, I);
    	z0 = z0(I);
        
	    init_a = zeros(1, K);
    	init_theta = zeros(1, K);
    	init_c = zeros(1, 2*(J-1)*K);
    	kk = 1;
    	for k = 1:K
        	init_a(k) = abs(z0(kk));
    	    init_theta(k) = abs(angle(z0(kk)));
    	    for j = 1:J-1
        	    init_c(2*(J-1)*(k-1)+2*j-1) = real(V0(j+1,kk) / V0(1,kk));
            	if imag(z0(kk)) < 0
                	init_c(2*(J-1)*(k-1)+2*j) = imag(V0(j+1,kk) / V0(1,kk));
	            else
    	            init_c(2*(J-1)*(k-1)+2*j) = -imag(V0(j+1,kk) / V0(1,kk));
        	    end
	        end
    	    if imag(z0(kk)) == 0, kk = kk + 1; else, kk = kk + 2; end
        end

        m = floor(T/2) + 1;
        freq0 = 2*pi*(0:m-1)'/T;
        P = zeros(J, J, length(freq0), K+1);
	    for k = 1:K
    	    a = init_a(k);
        	theta = init_theta(k);
            H = [1 0; reshape(init_c(2*(J-1)*(k-1)+1:2*(J-1)*k), J-1, 2)];
        	for i = 1:length(freq0)
                A = abs(1 - a*exp(1i*(theta - freq0(i))))^-2;
                B = abs(1 - a*exp(1i*(theta + freq0(i))))^-2;
            	P(:,:,i,k) = H * [A+B, 1i*(A-B); -1i*(A-B), A+B] * H' / (4*pi);
    	    end
        end
     	for i = 1:length(freq0)
            P(:,:,i,K+1) = eye(J) / (2*pi);
        end
        Y_fft = fft(Y, [], 2);
	    p_full = zeros(J, J, length(freq0));
    	for i = 1:length(freq0)
    	    p_full(:,:,i) = (Y_fft(:, i) * Y_fft(:, i)') / T;
        end
        options = optimoptions('fminunc','Algorithm','quasi-newton');
        beta = fminunc(@(b)whittle_regression_multi(b,P,p_full),zeros(size(P,4),1),options);
        weight = exp(beta);
        init_sigma2 = weight(1:K)';
        init_tau2 = weight(K+1);
        
        param0 = [atanh(2*init_a-1) zeros(1,K) log(init_sigma2/init_tau2) init_c]';

        switch lower(obj_func)
            case 'kalman'
                obj = @(param)kalman_multi_ll_wrapper(param,init_theta,Y);

            case 'whittle'
                freq_sub = freq0(1:ss_interval:end);
                p_sub = p_full(:,:,1:ss_interval:end);
                obj = @(param)whittle_multi_ll_wrapper(param, init_theta, freq_sub, p_sub);

            otherwise
                error('osc_decomp_multi: Unknown obj_func specified.')
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
                    p_batch = p_full(:,:,idx);
                    
                    [mll_t, ~, grad] = whittle_multi_ll(param, init_theta, freq_batch, p_batch);
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
                
            otherwise
                error('osc_decomp_multi: Unsupported optimization algorithm.');
        end
        
        param = param(:)';
    	[mll,osc_tau2] = kalman_multi_ll(param,init_theta,Y);
    	osc_AIC(K) = 2*mll + 2*((2*J+1)*K+1);
        param(K+1:2*K) = init_theta + tanh(param(K+1:2*K))*pi;
	    [~,I] = sort(abs(param(K+1:2*K)));
    	osc_a = (tanh(param(I))+1)/2;
    	osc_f = abs(param(K+I)) * fs / (2*pi);
    	osc_sigma2 = exp(param(2*K+I)) * osc_tau2;
        
        tmp = zeros(1, 2*(J-1)*K);
        for k = 1:K
            tmp((k-1)*2*(J-1)+1 : k*2*(J-1)) = I(k) - 1;
        end
        osc_c = param(3*K + tmp*2*(J-1) + repmat(1:2*(J-1), 1, K));
		
        osc_param(K, 1:(2*J+1)*K+1) = [osc_a osc_f osc_sigma2 osc_c osc_tau2];
        
        [osc_mean(1:2*K,:,K), osc_cov(1:2*K,1:2*K,:,K)] = osc_smooth(Y, fs, osc_param(K, 1:(2*J+1)*K+1));
        
	    for k = 1:K
        	osc_phase(k,:,K) = atan2(osc_mean(2*k,:,K), osc_mean(2*k-1,:,K));
    	end
	end
end

function [obj, grad] = whittle_regression_multi(b, X, P)
    b = b(:);
    B = length(b);
    [J, ~, M, ~] = size(X);
    
    exp_b = exp(b);    
    S = zeros(J, J, M);
    for m = 1:M
        S_m = zeros(J, J);
        for idx_b = 1:B
            S_m = S_m + X(:, :, m, idx_b) * exp_b(idx_b);
        end
        S_m = (S_m + S_m') / 2;
        if rcond(S_m) < 1e-12
            S_m = S_m + 1e-6 * eye(J);
        end
        S(:, :, m) = S_m;
    end
    
    obj = 0;
    for m = 1:M
        S_m = S(:, :, m);
        P_m = P(:, :, m);
        P_m = (P_m + P_m') / 2;        
        obj = obj + log(real(det(S_m))) + real(trace(S_m \ P_m));
    end
    
    if nargout > 1
        grad = zeros(B, 1);
        for idx_b = 1:B
            grad_b = 0;
            for m = 1:M
                S_m = S(:, :, m);
                P_m = P(:, :, m);
                P_m = (P_m + P_m') / 2;
                
                Inv_S = inv(S_m);
                W_m = Inv_S - (Inv_S * P_m * Inv_S);
                
                X_mb = X(:, :, m, idx_b);
                grad_b = grad_b + real(trace(W_m * X_mb));
            end
            grad(idx_b) = grad_b * exp_b(idx_b);
        end
    end
end

function [mll, grad] = kalman_multi_ll_wrapper(param, init_theta, Y)
    if nargout > 1
        [mll, ~, grad] = kalman_multi_ll(param, init_theta, Y);
    else
        mll = kalman_multi_ll(param, init_theta, Y);
    end
end

function [mll, grad] = whittle_multi_ll_wrapper(param, init_theta, freq, P)
    if nargout > 1
        [mll, ~, grad] = whittle_multi_ll(param, init_theta, freq, P);
    else
        mll = whittle_multi_ll(param, init_theta, freq, P);
    end
end
