function [mll, tau2hat, grad] = whittle_multi_ll(param, init_theta, freq, P)
    param = param(:);
    J = size(P, 1); 
    M = length(freq); 
    K = length(param) / (3 + 2*(J-1));
    
    u = param(1:K)';
    v = param(K+1:2*K)';
    w = param(2*K+1:3*K)';
    c_param = param(3*K+1:end); 
    
    a = (tanh(u) + 1) / 2;
    theta = init_theta(:)' + tanh(v) * pi;
    sigma2_ratio = exp(w);
    
    H = zeros(J, 2*K);
    H(1, 1:2:2*K) = 1;
    kk = 1;
    for k = 1:K
        for j = 2:J
            H(j, 2*k-1:2*k) = c_param(kk:kk+1)';
            kk = kk + 2;
        end
    end
    
    R_core = eye(J);
    I_2 = eye(2);
    
    S_y_ratio_inv = zeros(J, J, M);
    InvM_all = zeros(2, 2, K, M);
    
    log_det_sum = 0;
    tr_inv_P_sum = 0;
    
    for m_idx = 1:M
        omega = freq(m_idx);
        e_jomega = exp(-1i * omega);
        S_x_total = zeros(2*K, 2*K);
        
        for k = 1:K
            idx = 2*k-1:2*k;
            cos_t = cos(theta(k));
            sin_t = sin(theta(k));
            F_k = a(k) * [cos_t -sin_t; sin_t cos_t];
            
            InvTrans = I_2 - F_k * e_jomega;
            det_inv = InvTrans(1,1)*InvTrans(2,2) - InvTrans(1,2)*InvTrans(2,1);
            if abs(det_inv) < 1e-12, det_inv = det_inv + 1e-9; end
            InvM = [InvTrans(2,2), -InvTrans(1,2); -InvTrans(2,1), InvTrans(1,1)] / det_inv;
            
            InvM_all(:,:,k,m_idx) = InvM;
            S_x_total(idx, idx) = sigma2_ratio(k) * (InvM * InvM');
        end
        
        S_y_m = H * S_x_total * H' + R_core;
        S_y_m = (S_y_m + S_y_m') / 2;
        if rcond(S_y_m) < 1e-12, S_y_m = S_y_m + 1e-6 * eye(J); end
        
        S_y_ratio_inv(:,:,m_idx) = inv(S_y_m);
        
        P_m = P(:,:,m_idx);
        P_m = (P_m + P_m') / 2;
        
        tr_inv_P_sum = tr_inv_P_sum + real(trace(S_y_ratio_inv(:,:,m_idx) * P_m));
        log_det_sum = log_det_sum + log(real(det(S_y_m)));
    end
    
    tau2hat = tr_inv_P_sum / (J * M);
    if tau2hat < 1e-12 || isnan(tau2hat)
        tau2hat = 1e-9;
    end
    
    ll = -M * J * log(tau2hat) - log_det_sum - M * J;
    mll = -ll;
    
    if nargout < 3
        return;
    end
    
    grad_u = zeros(K, 1);
    grad_v = zeros(K, 1);
    grad_w = zeros(K, 1);
    grad_c = zeros(length(c_param), 1);
    
    W = zeros(J, J, M);
    for m_idx = 1:M
        Inv_S = S_y_ratio_inv(:,:,m_idx);
        W(:,:,m_idx) = (Inv_S * P(:,:,m_idx) * Inv_S) / tau2hat - Inv_S;
    end
    
    for k = 1:K
        idx = 2*k-1:2*k;
        H_k = H(:, idx); 
        
        da_du = 2 * a(k) * (1 - a(k));
        dtheta_dv = pi * (1 - tanh(v(k))^2);
        dsigma2_dw = sigma2_ratio(k);
        
        dl_da = 0;
        dl_dtheta = 0;
        dl_dsigma2 = 0;
        
        for m_idx = 1:M
            omega = freq(m_idx);
            e_jomega = exp(-1i * omega);
            InvM = InvM_all(:,:,k,m_idx);
            W_m = W(:,:,m_idx);
            
            HwH = H_k' * W_m * H_k; 
            
            dSx_dsigma2 = InvM * InvM';
            dl_dsigma2 = dl_dsigma2 + real(trace(HwH * dSx_dsigma2));
            
            cos_t = cos(theta(k));
            sin_t = sin(theta(k));
            
            dF_da = [cos_t -sin_t; sin_t cos_t];
            dF_dtheta = a(k) * [-sin_t -cos_t; cos_t -sin_t];
            
            dInvM_da = e_jomega * (InvM * dF_da * InvM);
            dInvM_dtheta = e_jomega * (InvM * dF_dtheta * InvM);
            
            dSx_da = sigma2_ratio(k) * (dInvM_da * InvM' + InvM * dInvM_da');
            dSx_dtheta = sigma2_ratio(k) * (dInvM_dtheta * InvM' + InvM * dInvM_dtheta');
            
            dl_da = dl_da + real(trace(HwH * dSx_da));
            dl_dtheta = dl_dtheta + real(trace(HwH * dSx_dtheta));
        end
        
        grad_u(k) = dl_da * da_du;
        grad_v(k) = dl_dtheta * dtheta_dv;
        grad_w(k) = dl_dsigma2 * dsigma2_dw;
    end
    
    kk = 1;
    for k = 1:K
        idx = 2*k-1:2*k;
        dH_accum = zeros(J, 2);
        for m_idx = 1:M
            InvM = InvM_all(:,:,k,m_idx);
            S_x_k = sigma2_ratio(k) * (InvM * InvM');
            dH_accum = dH_accum + 2 * (W(:,:,m_idx) * H(:, idx) * S_x_k);
        end
        
        for j = 2:J
            grad_c(kk)   = real(dH_accum(j, 1)); 
            grad_c(kk+1) = real(dH_accum(j, 2)); 
            kk = kk + 2;
        end
    end
    
    grad = [grad_u; grad_v; grad_w; grad_c];
    grad = -grad;
    grad(isnan(grad)) = 0;
end
