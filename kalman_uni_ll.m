function [mll, tau2hat, grad] = kalman_uni_ll(param,init_theta,y)
    T = length(y);
    K = length(param) / 3;
    
    a = (tanh(param(1:K)) + 1) / 2;
    theta = init_theta + tanh(param(K+1:2*K)) * pi;
    sigma2 = exp(param(2*K+1:3*K));
    
    F = zeros(2*K, 2*K);
    Q = zeros(2*K, 2*K);
    for k = 1:K
        F(2*k-1:2*k, 2*k-1:2*k) = a(k) * [cos(theta(k)) -sin(theta(k)); sin(theta(k)) cos(theta(k))];
        Q(2*k-1:2*k, 2*k-1:2*k) = sigma2(k) * eye(2);
    end
    
    H = zeros(1, 2*K);
    H(1:2:2*K) = 1;
    R = 1; 
    
    x_pred1 = zeros(2*K, T);
    x_filt  = zeros(2*K, T);
    V_pred1 = zeros(2*K, 2*K, T);
    V_filt  = zeros(2*K, 2*K, T);
    
    x_pred1(:, 1) = zeros(2*K, 1);
    for k = 1:K
        V_pred1(2*k-1:2*k, 2*k-1:2*k, 1) = sigma2(k) / (1 - a(k)^2) * eye(2);
    end
    
    err_var_vec = zeros(1, T); 
    err_vec = zeros(1, T);
    I_state = eye(2*K);
    
    for t = 1:T
        V_p = V_pred1(:,:,t);
        V_p = (V_p + V_p') / 2;
        
        err_var = H * V_p * H' + R;
        if err_var < 1e-12 || isnan(err_var)
            err_var = err_var + 1e-9;
        end
        err_var_vec(t) = err_var;
        
        err_vec(t) = y(t) - H * x_pred1(:,t);
        Kg = (V_p * H') / err_var; 
        
        x_filt(:,t) = x_pred1(:,t) + Kg * err_vec(t);
        
        ImKH = I_state - Kg * H;
        V_f = ImKH * V_p * ImKH' + Kg * R * Kg';
        V_f = (V_f + V_f') / 2;
        V_filt(:,:,t) = V_f;
        
        if t < T
            x_pred1(:,t+1) = F * x_filt(:,t);
            V_p_next = F * V_f * F' + Q;
            V_pred1(:,:,t+1) = (V_p_next + V_p_next') / 2;
        end
    end
    
    tau2hat = sum((err_vec.^2) ./ err_var_vec) / T;
    ll = - (T / 2) * log(tau2hat) - (T / 2) - (T / 2) * log(2*pi) - 0.5 * sum(log(err_var_vec));
    mll = -ll;
    
    if nargout < 3
        return;
    end
    
    grad_F = zeros(2*K, 2*K, 3*K);
    grad_Q = zeros(2*K, 2*K, 3*K);
    for k = 1:K
        grad_F(2*k-1:2*k, 2*k-1:2*k, k) = (1 / (2 * cosh(a(k))^2)) * [cos(theta(k)) -sin(theta(k)); sin(theta(k)) cos(theta(k))];
        grad_F(2*k-1:2*k, 2*k-1:2*k, K+k) = (pi * a(k) / cosh(param(K+k))^2) * [-sin(theta(k)) -cos(theta(k)); cos(theta(k)) -sin(theta(k))];
        grad_Q(2*k-1:2*k, 2*k-1:2*k, 2*K+k) = sigma2(k) * eye(2);
    end
    
    grad_x_pred1 = zeros(2*K, 3*K);
    grad_V_pred1 = zeros(2*K, 2*K, 3*K);
    for k = 1:K
        grad_V_pred1(2*k-1:2*k, 2*k-1:2*k, k) = (a(k) * sigma2(k)) / (cosh(a(k))^2 * (1 - a(k)^2)^2) * eye(2);
        grad_V_pred1(2*k-1:2*k, 2*k-1:2*k, 2*K+k) = sigma2(k) / (1 - a(k)^2) * eye(2);
    end
    
    grad_x_filt = zeros(2*K, 3*K);
    grad_V_filt = zeros(2*K, 2*K, 3*K);
    grad_Rhat   = zeros(3*K, 1);
    grad        = zeros(3*K, 1);
    
    for t = 1:T
        err = err_vec(t);
        err_var = err_var_vec(t);
        Kg = (V_pred1(:,:,t) * H') / err_var;
        
        grad_err = zeros(3*K, 1);
        grad_err_var = zeros(3*K, 1);
        for i = 1:3*K
            grad_err(i) = -H * grad_x_pred1(:,i);
            grad_err_var(i) = H * grad_V_pred1(:,:,i) * H';
        end
        
        for i = 1:3*K
            grad(i) = grad(i) - (grad_err_var(i) / err_var) / 2;
            grad_Rhat(i) = grad_Rhat(i) - 2 * err * (-grad_err(i)) / (err_var * T) ...
                                        - (err^2 / err_var^2) * grad_err_var(i) / T;
        end
        
        grad_Kg = zeros(2*K, 3*K);
        for i = 1:3*K
            grad_Kg(:,i) = (grad_V_pred1(:,:,i) * H') / err_var - Kg * grad_err_var(i) / err_var;
        end
        
        for i = 1:3*K
            grad_x_filt(:,i) = grad_x_pred1(:,i) + Kg * grad_err(i) + grad_Kg(:,i) * err;
            
            ImKH = I_state - Kg * H;
            d_ImKH = -grad_Kg(:,i) * H;
            grad_V_filt(:,:,i) = d_ImKH * V_pred1(:,:,t) * ImKH' + ImKH * grad_V_pred1(:,:,i) * ImKH' ...
                                 + ImKH * V_pred1(:,:,t) * d_ImKH' + grad_Kg(:,i) * R * Kg' + Kg * R * grad_Kg(:,i)';
            grad_V_filt(:,:,i) = (grad_V_filt(:,:,i) + grad_V_filt(:,:,i)') / 2;
        end
        
        if t < T
            for i = 1:3*K
                grad_x_pred1(:,i) = F * grad_x_filt(:,i) + grad_F(:,:,i) * x_filt(:,t);
                
                grad_V_p_next = F * grad_V_filt(:,:,i) * F' + grad_F(:,:,i) * V_filt(:,:,t) * F' ...
                                + F * V_filt(:,:,t) * grad_F(:,:,i)' + grad_Q(:,:,i);
                grad_V_pred1(:,:,i) = (grad_V_p_next + grad_V_p_next') / 2;
            end
        end
    end
    
    grad = grad - (T / (2 * tau2hat)) * grad_Rhat;
    grad = -grad;
end
