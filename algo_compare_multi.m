MAX_OSC = 5;
J = 2;
K = 3;
T = 10000;
fs = 200;
a = [0.9 0.9 0.9];
f = [20 50 75];
sigma2 = [.1 .1 .1];
c = [1.8 1.2 0 1 -0.1 0.4];
tau2 = 0.1;
x = zeros(2*K,T);
for k=1:K
    x(2*k-1:2*k,1) = sqrt(sigma2(k))/sqrt(1-a(k)^2)*randn(2,1);
end
for t=2:T
    for k=1:K
        x(2*k-1:2*k,t) = a(k)*[cos(2*pi*f(k)/fs) -sin(2*pi*f(k)/fs); sin(2*pi*f(k)/fs) cos(2*pi*f(k)/fs)]*x(2*k-1:2*k,t-1)+sqrt(sigma2(k))*randn(2,1);
    end
end
H = zeros(J,2*K);
H(1,1:2:2*K) = 1;
H(2,:) = c;
Y = H*x+sqrt(tau2)*randn(J,T);

tic
[osc_param_kalman_qn0,osc_AIC_kalman_qn0] = osc_decomp_multi(Y,fs,MAX_OSC,'quasi-newton_grad0','kalman');
time_kalman_qn0 = toc
osc_AIC_kalman_qn0
tic
[osc_param_kalman_qn1,osc_AIC_kalman_qn1] = osc_decomp_multi(Y,fs,MAX_OSC,'quasi-newton_grad1','kalman');
time_kalman_qn1 = toc
osc_AIC_kalman_qn1
tic
[osc_param_whittle_qn0,osc_AIC_whittle_qn0] = osc_decomp_multi(Y,fs,MAX_OSC,'quasi-newton_grad0','whittle');
time_whittle_qn0 = toc
osc_AIC_whittle_qn0
tic
[osc_param_whittle_qn1,osc_AIC_whittle_qn1] = osc_decomp_multi(Y,fs,MAX_OSC,'quasi-newton_grad1','whittle');
time_whittle_qn1 = toc
osc_AIC_whittle_qn1
tic
[osc_param_whittle_qn0_100,osc_AIC_whittle_qn0_100] = osc_decomp_multi(Y,fs,MAX_OSC,'quasi-newton_grad0','whittle',100);
time_whittle_qn0_100 = toc
osc_AIC_whittle_qn0_100
tic
[osc_param_whittle_qn1_100,osc_AIC_whittle_qn1_100] = osc_decomp_multi(Y,fs,MAX_OSC,'quasi-newton_grad1','whittle',100);
time_whittle_qn1_100 = toc
osc_AIC_whittle_qn1_100
tic
[osc_param_whittle_adam1,osc_AIC_whittle_adam1] = osc_decomp_multi(Y,fs,MAX_OSC,'adam','whittle',1);
time_whittle_adam1 = toc
osc_AIC_whittle_adam1
tic
[osc_param_whittle_adam10,osc_AIC_whittle_adam10] = osc_decomp_multi(Y,fs,MAX_OSC,'adam','whittle',10);
time_whittle_adam10 = toc
osc_AIC_whittle_adam10
tic
[osc_param_whittle_adam100,osc_AIC_whittle_adam100] = osc_decomp_multi(Y,fs,MAX_OSC,'adam','whittle',100);
time_whittle_adam100 = toc
osc_AIC_whittle_adam100
save('algo_compare_multi')
