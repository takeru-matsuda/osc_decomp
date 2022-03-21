% Kalman smoother for oscillator model (both univariate and multivariate)
function [x_smooth,V_smooth] = osc_smooth(Y,fs,param)
    J = size(Y,1);
    T = size(Y,2);
    K = (length(param)-1)/(3+2*(J-1));
    a = param(1:K);
    theta  = param(K+1:2*K)/fs*2*pi;
    sigma2 = param(2*K+1:3*K);
    c = param(3*K+1:end-1);
    tau2 = param(end);
    
    F = zeros(2*K,2*K);
    Q = zeros(2*K,2*K);
    for k=1:K
        F(2*k-1:2*k,2*k-1:2*k) = a(k)*[cos(theta(k)) -sin(theta(k)); sin(theta(k)) cos(theta(k))];
        Q(2*k-1:2*k,2*k-1:2*k) = sigma2(k)*eye(2);
    end
    H = zeros(J,2*K);
    H(1,1:2:2*K) = 1;
    kk = 1;
    for k=1:K
        for j=2:J
            H(j,2*k-1:2*k) = c(kk:kk+1);
            kk = kk+2;
        end
    end
    R = tau2*eye(J);
    
    x_pred1 = zeros(2*K,T);
    x_filt = zeros(2*K,T);
    x_smooth = zeros(2*K,T);
    V_pred1 = zeros(2*K,2*K,T);
    V_filt = zeros(2*K,2*K,T);
    x_pred1(:,1) = zeros(2*K,1);
    for k=1:K
        V_pred1(2*k-1:2*k,2*k-1:2*k,1) = sigma2(k)/(1-a(k)^2)*eye(2);
    end
    for t=1:T
        x_filt(:,t) = x_pred1(:,t) + V_pred1(:,:,t)*H'*((H*V_pred1(:,:,t)*H'+R)\(Y(:,t)-H*x_pred1(:,t)));
        V_filt(:,:,t) = V_pred1(:,:,t) - V_pred1(:,:,t)*H'*((H*V_pred1(:,:,t)*H'+R)\H)*V_pred1(:,:,t);
        if t == T
            break
        end
        x_pred1(:,t+1) = F*x_filt(:,t);
        V_pred1(:,:,t+1) = F*V_filt(:,:,t)*F'+Q;
    end
    x_smooth(:,T) = x_filt(:,T);
    V_smooth(:,:,T) = V_filt(:,:,T);
    for t=T-1:-1:1
        x_smooth(:,t) = x_filt(:,t) + V_filt(:,:,t)*F'*(V_pred1(:,:,t+1)\(x_smooth(:,t+1)-x_pred1(:,t+1)));
        V_smooth(:,:,t) = V_filt(:,:,t) + V_filt(:,:,t)*F'*(V_pred1(:,:,t+1)\(V_smooth(:,:,t+1)-V_pred1(:,:,t+1)))*(V_pred1(:,:,t+1)\F)*V_filt(:,:,t);
    end
end


