function [x_smooth,V_smooth] = osc_smooth(Y,fs,param)
    param = param(:)';
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
    V_smooth = zeros(2*K,2*K,T);
    
    x_pred1(:,1) = zeros(2*K,1);
    for k=1:K
        V_pred1(2*k-1:2*k,2*k-1:2*k,1) = sigma2(k)/(1-a(k)^2)*eye(2);
    end
    for t=1:T
        V_p_t = V_pred1(:,:,t);
        V_p_t = (V_p_t + V_p_t') / 2;
        InnCov = H*V_p_t*H' + R;
        InnCov = (InnCov + InnCov') / 2;
        if rcond(InnCov) < 1e-12
            InnCov = InnCov + 1e-9 * eye(size(InnCov));
        end
        K_g = (InnCov \ (H * V_p_t))';
        x_filt(:,t) = x_pred1(:,t) + K_g * (Y(:,t) - H*x_pred1(:,t));
        ImKH = eye(2*K) - K_g * H;
        V_f_t = ImKH * V_p_t * ImKH' + K_g * R * K_g';
        V_f_t = (V_f_t + V_f_t') / 2;
        V_filt(:,:,t) = V_f_t;        
        if t == T
            break
        end
        x_pred1(:,t+1) = F*x_filt(:,t);
        V_pred1(:,:,t+1) = F*V_filt(:,:,t)*F'+Q;
    end

    x_smooth(:,T) = x_filt(:,T);
    V_smooth(:,:,T) = V_filt(:,:,T);
    for t=T-1:-1:1
        V_p = V_pred1(:,:,t+1);
        V_p = (V_p + V_p') / 2; 
        if rcond(V_p) < 1e-12 || isnan(rcond(V_p))
            V_p = V_p + 1e-9 * eye(size(V_p));
        end
        J_t = (V_p \ (F * V_filt(:,:,t)'))';        
        x_smooth(:,t) = x_filt(:,t) + J_t * (x_smooth(:,t+1) - x_pred1(:,t+1));
        V_smooth(:,:,t) = V_filt(:,:,t) + J_t * (V_smooth(:,:,t+1) - V_p) * J_t';        
        V_smooth(:,:,t) = (V_smooth(:,:,t) + V_smooth(:,:,t)') / 2;    
    end
end
