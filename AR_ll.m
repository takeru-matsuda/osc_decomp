function [mll,Ehat] = AR_ll(y,p)
    T = length(y);
    ARdeg = length(p);
    c = (exp(p)-1)./(exp(p)+1);
    a = zeros(ARdeg,ARdeg);
    for m=1:ARdeg
        a(m,m) = c(m);
    end
    for m=2:ARdeg
        for i=1:m-1
            a(i,m) = a(i,m-1)-c(m)*a(m-i,m-1);
        end
    end
    K = zeros(ARdeg+1,ARdeg+1);
    A = [1 -a(:,ARdeg)'];
    for i=1:ARdeg+1
        K(i,i:-1:2) = K(i,i:-1:2)+A(1:i-1);
        K(i,1:ARdeg-i+2) = K(i,1:ARdeg-i+2)+A(i:ARdeg+1);
    end
    C = K\[1; zeros(ARdeg,1)];
    F = zeros(ARdeg,ARdeg);
    F(:,1) = a(:,ARdeg);
    F(1:ARdeg-1,2:ARdeg) = eye(ARdeg-1);
    G = [1; zeros(ARdeg-1,1)];
    Q = G*G';
    H = [1 zeros(1,ARdeg-1)];
    R = 0;
    x_pred1 = zeros(ARdeg,T);
    x_filt = zeros(ARdeg,T);
    V_pred1 = zeros(ARdeg,ARdeg,T);
    V_filt = zeros(ARdeg,ARdeg,T);
    x_pred1(:,1) = zeros(ARdeg,1);
    V_pred1(1,1,1) = C(1);
    for i=2:ARdeg
        V_pred1(i,1,1) = C(2:ARdeg-i+2)'*a(i:ARdeg,ARdeg);
        V_pred1(1,i,1) = V_pred1(i,1,1);
    end
    for i=2:ARdeg
        for j=i:ARdeg
            for p=i:ARdeg
                for q=j:ARdeg
                    V_pred1(i,j,1) = V_pred1(i,j,1)+a(p,ARdeg)*a(q,ARdeg)*C(abs(q-j-p+i)+1);
                end
            end
            V_pred1(j,i,1) = V_pred1(i,j,1);
        end
    end
    for t=1:T-1
        x_filt(:,t) = x_pred1(:,t) + V_pred1(:,:,t)*H'*((H*V_pred1(:,:,t)*H'+R)\(y(t)-H*x_pred1(:,t)));
        V_filt(:,:,t) = (eye(ARdeg)-V_pred1(:,:,t)*(H'*H)/(H*V_pred1(:,:,t)*H'+R))*V_pred1(:,:,t);
        x_pred1(:,t+1) = F*x_filt(:,t);
        V_pred1(:,:,t+1) = F*V_filt(:,:,t)*F'+Q;
    end
    Ehat = 0;
    for t=1:T
        Ehat = Ehat*(t-1)/t+(y(t)-H*x_pred1(:,t))^2/(H*V_pred1(:,:,t)*H')/t;
    end
    ll = -T*log(Ehat)/2-T/2-T/2*log(2*pi);
    for t=1:T
        ll = ll-log(H*V_pred1(:,:,t)*H')/2;
    end
    mll = -ll;
end

