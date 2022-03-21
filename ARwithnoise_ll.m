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
    for t=1:T-1
        x_filt(:,t) = x_pred1(:,t) + V_pred1(:,:,t)*H'*((H*V_pred1(:,:,t)*H'+R)\(y(t)-H*x_pred1(:,t)));
        V_filt(:,:,t) = V_pred1(:,:,t) - V_pred1(:,:,t)*H'*((H*V_pred1(:,:,t)*H'+R)\H)*V_pred1(:,:,t);
        x_pred1(:,t+1) = F*x_filt(:,t);
        V_pred1(:,:,t+1) = F*V_filt(:,:,t)*F'+Q;
    end
    Rhat = 0;
    for t=1:T
        Rhat = Rhat*(t-1)/t+(y(t)-H*x_pred1(:,t))^2/(H*V_pred1(:,:,t)*H'+R)/t;
    end
    ll = -T*log(Rhat)/2-T/2-T/2*log(2*pi);
    for t=1:T
        ll = ll-log(H*V_pred1(:,:,t)*H'+R)/2;
    end
    % for minimization
    mll = -ll;
end

