function [mll,Rhat] = ARwithnoise_ll(y,param)
    T = length(y);
    ARdeg = length(param)-1;
    A = [1 param(1:ARdeg)];
    E = exp(param(ARdeg+1));
    R = 1;
    x_pred1 = zeros(ARdeg,T);
    x_filt = zeros(ARdeg,T);
    V_pred1 = zeros(ARdeg,ARdeg,T);
    V_filt = zeros(ARdeg,ARdeg,T);
    F = [-A(2:ARdeg+1); eye(ARdeg-1) zeros(ARdeg-1,1)];
    Q = [E zeros(1,ARdeg-1); zeros(ARdeg-1,ARdeg)];
    H = [1 zeros(1,ARdeg-1)];
    x_pred1(:,1) = zeros(ARdeg,1);
    K = zeros(ARdeg+1,ARdeg+1);
    for i=1:ARdeg+1
        K(i,i:-1:2) = K(i,i:-1:2)+A(1:i-1);
        K(i,1:ARdeg-i+2) = K(i,1:ARdeg-i+2)+A(i:ARdeg+1);
    end
    c = K\[E; zeros(ARdeg,1)];
    V_pred1(:,:,1) = toeplitz(c(1:ARdeg));
    
    err_var_vec = zeros(1,T);
    err_vec = zeros(1,T);
    
    for t=1:T-1
        err_vec(t) = y(t) - H*x_pred1(:,t);
        err_var = H*V_pred1(:,:,t)*H' + R;
        err_var_vec(t) = err_var;
        
        Kg = (V_pred1(:,:,t)*H') / err_var;
        
        x_filt(:,t) = x_pred1(:,t) + Kg * err_vec(t);
        V_filt(:,:,t) = V_pred1(:,:,t) - Kg * H * V_pred1(:,:,t);
        x_pred1(:,t+1) = F*x_filt(:,t);
        V_pred1(:,:,t+1) = F*V_filt(:,:,t)*F'+Q;
    end
    
    err_vec(T) = y(T) - H*x_pred1(:,T);
    err_var_vec(T) = H*V_pred1(:,:,T)*H' + R;
    
    Rhat = sum((err_vec.^2) ./ err_var_vec) / T;
    ll = - (T/2)*log(Rhat) - T/2 - (T/2)*log(2*pi) - 0.5*sum(log(err_var_vec));
    mll = -ll;
end