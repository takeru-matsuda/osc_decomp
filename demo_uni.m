load('CanadianLynxData');
y = log(lynx);  % not necessary for neural data
y = y-mean(y);
fs = 1; % sampling frequency
MAX_OSC = 6;
[osc_param,osc_AIC,osc_mean,osc_cov,osc_phase] = osc_decomp(y,fs,MAX_OSC);
[minAIC,K] = min(osc_AIC);
osc_a = osc_param(K,1:K);
osc_f = osc_param(K,K+1:2*K);
osc_sigma2 = osc_param(K,2*K+1:3*K);
osc_tau2 = osc_param(K,3*K+1);
[hess,grad,mll] = kalman_ll_hess(y,fs,osc_param(K,1:3*K+1));
[cov_est,ci_info] = osc_hessian_covariance(hess);
if ~ci_info.valid
    fprintf('Parameter confidence intervals unavailable: %s.\n',ci_info.reason);
end
fprintf('The number of oscillators is K=%d.\n',K);
fprintf('The periods of K oscillators are:\n');
for k=1:K
    if osc_f(k)==0
        fprintf(' Component %d: zero frequency; period undefined.\n',k);
    elseif ci_info.valid
        fprintf(' %.2f (95%% CI: [%.2f %.2f]) years\n',1./osc_f(k),max(0,1./osc_f(k)-1.96*sqrt(cov_est(K+k,K+k))/osc_f(k)^2),1./osc_f(k)+1.96*sqrt(cov_est(K+k,K+k))/osc_f(k)^2);
    else
        fprintf(' %.2f years (95%% CI unavailable).\n',1./osc_f(k));
    end
end
osc_plot(osc_mean,osc_cov,fs,K)
osc_phase_plot(osc_phase,osc_mean,osc_cov,fs,K)
osc_spectrum_plot(y,fs,osc_a,osc_f,osc_sigma2,osc_tau2)

