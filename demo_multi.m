load('NorthSouthSunspotData');
Y = log(dat+1);
Y = Y-mean(Y,2)*ones(1,size(Y,2));
fs = 12;
J = size(Y,1);
T = size(Y,2);
MAX_OSC = 6;
[osc_param,osc_AIC,osc_mean,osc_cov,osc_phase] = osc_decomp(Y,fs,MAX_OSC);
[minAIC,K] = min(osc_AIC);
osc_a = osc_param(K,1:K);
osc_f = osc_param(K,K+1:2*K);
osc_sigma2 = osc_param(K,2*K+1:3*K);
osc_c = osc_param(K,3*K+1:(2*J+1)*K);
osc_tau2 = osc_param(K,(2*J+1)*K+1);
hess = osc_ll_hess(Y,fs,osc_param(K,1:(2*J+1)*K+1));
cov_est = inv(hess);
fprintf('The number of oscillators is K=%d.\n',K);
fprintf('The periods of K oscillators are:\n');
for k=1:K
    fprintf(' %.2f (95%% CI: [%.2f %.2f]) years\n',1./osc_f(k),1./(osc_f(k)+1.96*sqrt(cov_est(K+k,K+k))),1./(osc_f(k)-1.96*sqrt(cov_est(K+k,K+k))));
end
fprintf('The phase differences for K oscillators correspond to:\n');
for k=1:K
    phase_diff = atan2(osc_c(2*k),osc_c(2*k-1));
    tmp = [-osc_c(2*k); osc_c(2*k-1)]/(osc_c(2*k-1)^2+osc_c(2*k)^2);
    phase_var_est = tmp'*cov_est(3*K+2*k-1:3*K+2*k,3*K+2*k-1:3*K+2*k)*tmp;
    fprintf(' %.2f (95%% CI: [%.2f %.2f]) years\n',phase_diff/2/pi/osc_f(k),(phase_diff-1.96*sqrt(phase_var_est))/2/pi/osc_f(k),(phase_diff+1.96*sqrt(phase_var_est))/2/pi/osc_f(k));
end
osc_plot(osc_mean,osc_cov,fs,K);
osc_phase_plot(osc_phase,osc_mean,osc_cov,fs,K);
osc_spectrum_plot(Y,fs,osc_a,osc_f,osc_sigma2,osc_tau2,osc_c);
