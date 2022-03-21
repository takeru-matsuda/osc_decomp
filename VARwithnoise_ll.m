% minus log-likelihood of stationary VAR model with observation noise
function mll = VARwithnoise_ll(Y,A,E,R)
    EPS = 10^-6;
    J = size(Y,1);
    T = size(Y,2);
    ARdeg = size(A,2)/J;
    x_pred1 = zeros(J*ARdeg,T);
    x_filt = zeros(J*ARdeg,T);
    V_pred1 = zeros(J*ARdeg,J*ARdeg,T);
    V_filt = zeros(J*ARdeg,J*ARdeg,T);
    F = [A(:,1:J*ARdeg); eye(J*(ARdeg-1)) zeros(J*(ARdeg-1),J)];
    Q = [E zeros(J,J*(ARdeg-1)); zeros(J*(ARdeg-1),J*ARdeg)];
    H = [eye(J) zeros(J,J*(ARdeg-1))];
    x_pred1(:,1) = zeros(J*ARdeg,1);
    V_pred1(:,:,1) = F*Q*F'+Q;
    prev = Q;
    cnt = 0;
    while norm(V_pred1(:,:,1)-prev,'fro')/norm(prev,'fro')>EPS
        prev = V_pred1(:,:,1);
        V_pred1(:,:,1) = F*V_pred1(:,:,1)*F'+Q;
        cnt = cnt+1;
    end
    for t=1:T-1
        x_filt(:,t) = x_pred1(:,t) + V_pred1(:,:,t)*H'*((H*V_pred1(:,:,t)*H'+R)\(Y(:,t)-H*x_pred1(:,t)));
        V_filt(:,:,t) = V_pred1(:,:,t) - V_pred1(:,:,t)*H'*((H*V_pred1(:,:,t)*H'+R)\H)*V_pred1(:,:,t);
        x_pred1(:,t+1) = F*x_filt(:,t);
        V_pred1(:,:,t+1) = F*V_filt(:,:,t)*F'+Q;
    end
    ll = -J*T/2*log(2*pi);
    for t=1:T
        ll = ll-log(det(H*V_pred1(:,:,t)*H'+R))/2-(Y(:,t)-H*x_pred1(:,t))'*((H*V_pred1(:,:,t)*H'+R)\(Y(:,t)-H*x_pred1(:,t)))/2;
    end
    % for minimization
    mll = -ll;
end

