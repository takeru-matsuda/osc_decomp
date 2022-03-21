T = 600;
K = 3;
fs = 200;
a = [0.9 0.9 0.9];
f = [20 50 75];
sigma2 = [.1 .1 .1];
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
H = zeros(1,2*K);
H(1:2:2*K) = 1;
y = H*x+sqrt(tau2)*randn(1,T);

MAX_OSC = 5;
[osc_param,osc_AIC,osc_mean,osc_cov,osc_phase] = osc_decomp(y,fs,MAX_OSC);
[minAIC,K] = min(osc_AIC);
osc_a = osc_param(K,1:K);
osc_f = osc_param(K,K+1:2*K);
osc_sigma2 = osc_param(K,2*K+1:3*K);
osc_tau2 = osc_param(K,3*K+1);
