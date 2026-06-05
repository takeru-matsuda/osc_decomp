function mll = VARwithnoise_ll(Y,A,E,R)
    EPS = 10^-6;
    J = size(Y,1);
    T = size(Y,2);
    ARdeg = size(A,2)/J;
    
    x_pred1 = zeros(J*ARdeg,T);
    x_filt = zeros(J*ARdeg,T);
    V_pred1 = zeros(J*ARdeg,J*ARdeg,T);
    V_filt = zeros(J*ARdeg,J*ARdeg,T);
    
    F = [A; eye(J*(ARdeg-1)) zeros(J*(ARdeg-1),J)];
    Q = [E zeros(J,J*(ARdeg-1)); zeros(J*(ARdeg-1),J*ARdeg)];
    H = [eye(J) zeros(J,J*(ARdeg-1))];
    
    V_p = F*Q*F'+Q;
    prev = Q;
    while norm(V_p - prev,'fro')/norm(prev,'fro') > EPS
        prev = V_p;
        V_p = F*V_p*F'+Q;
    end
    V_pred1(:,:,1) = V_p;
    
    ll_sum = 0;
    
    for t=1:T-1
        V_t = V_pred1(:,:,t);
        x_t = x_pred1(:,t);
        Y_err = Y(:,t) - H*x_t;
        S_t = H*V_t*H' + R;
        
        K_gain = V_t * H' / S_t;
        
        x_f = x_t + K_gain * Y_err;
        V_f = V_t - K_gain * H * V_t;
        
        x_filt(:,t) = x_f;
        V_filt(:,:,t) = V_f;
        
        x_pred1(:,t+1) = F*x_f;
        V_pred1(:,:,t+1) = F*V_f*F' + Q;
        
        ll_sum = ll_sum - log(det(S_t))/2 - (Y_err' / S_t * Y_err)/2;
    end
    
    V_T = V_pred1(:,:,T);
    Y_err_T = Y(:,T) - H*x_pred1(:,T);
    S_T = H*V_T*H' + R;
    ll_sum = ll_sum - log(det(S_T))/2 - (Y_err_T' / S_T * Y_err_T)/2;
    
    ll = -J*T/2*log(2*pi) + ll_sum;
    mll = -ll;
end
