function [hess,grad,mll] = osc_ll_hess(Y,fs,param)
    J = size(Y,1);
    T = size(Y,2);
    K = (length(param)-1)/(3+2*(J-1));
    a = param(1:K);
    f = param(K+1:2*K);
    sigma2 = param(2*K+1:3*K);
    c = param(3*K+1:end-1);
    tau2 = param(end);
    theta = 2*pi*f/fs;
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
        
    grad_F = zeros(2*K,2*K,length(param));
    grad_Q = zeros(2*K,2*K,length(param));
    grad_H = zeros(J,2*K,length(param));
    hess_F = zeros(2*K,2*K,length(param),length(param));
    kk = 1;
    for k=1:K
        grad_F(2*k-1:2*k,2*k-1:2*k,k) = [cos(theta(k)) -sin(theta(k)); sin(theta(k)) cos(theta(k))];
        grad_F(2*k-1:2*k,2*k-1:2*k,K+k) = a(k)*[-sin(theta(k)) -cos(theta(k)); cos(theta(k)) -sin(theta(k))];
        grad_Q(2*k-1:2*k,2*k-1:2*k,2*K+k) = eye(2);
        for j=2:J
            grad_H(j,2*k-1,3*K+kk) = 1;
            grad_H(j,2*k,3*K+kk+1) = 1;
            kk = kk+2;
        end
        hess_F(2*k-1:2*k,2*k-1:2*k,k,K+k) = [-sin(theta(k)) -cos(theta(k)); cos(theta(k)) -sin(theta(k))];
        hess_F(2*k-1:2*k,2*k-1:2*k,K+k,k) = [-sin(theta(k)) -cos(theta(k)); cos(theta(k)) -sin(theta(k))];
        hess_F(2*k-1:2*k,2*k-1:2*k,K+k,K+k) = a(k)*[-cos(theta(k)) sin(theta(k)); -sin(theta(k)) -cos(theta(k))];
    end
    grad_R = zeros(J,J,length(param));
    grad_R(:,:,end) = eye(J);

    x_pred1 = zeros(2*K,1);
    V_pred1 = zeros(2*K,2*K);
    x_filt = zeros(2*K,1);
    V_filt = zeros(2*K,2*K);
    grad_x_pred1 = zeros(2*K,length(param));
    grad_V_pred1 = zeros(2*K,2*K,length(param));
    grad_x_filt = zeros(2*K,length(param));
    grad_V_filt = zeros(2*K,2*K,length(param));
    hess_x_pred1 = zeros(2*K,length(param),length(param));
    hess_V_pred1 = zeros(2*K,2*K,length(param),length(param));
    hess_x_filt = zeros(2*K,length(param),length(param));
    hess_V_filt = zeros(2*K,2*K,length(param),length(param));
    for k=1:K
        V_pred1(2*k-1:2*k,2*k-1:2*k) = sigma2(k)/(1-a(k)^2)*eye(2);
        grad_V_pred1(2*k-1:2*k,2*k-1:2*k,k) = 2*a(k)*sigma2(k)/(1-a(k)^2)^2*eye(2);
        grad_V_pred1(2*k-1:2*k,2*k-1:2*k,2*K+k) = 1/(1-a(k)^2)*eye(2);
        hess_V_pred1(2*k-1:2*k,2*k-1:2*k,k,k) = 2*(3*a(k)^2+1)*sigma2(k)/(1-a(k)^2)^3*eye(2);
        hess_V_pred1(2*k-1:2*k,2*k-1:2*k,k,2*K+k) = 2*a(k)/(1-a(k)^2)^2*eye(2);
        hess_V_pred1(2*k-1:2*k,2*k-1:2*k,2*K+k,k) = 2*a(k)/(1-a(k)^2)^2*eye(2);
    end
    ll = -J*T/2*log(2*pi);
    grad = zeros(length(param),1);
    hess = zeros(length(param),length(param));
    for t=1:T
        err = Y(:,t)-H*x_pred1;
        err_cov = H*V_pred1*H'+R;
        inv_err_cov = inv(err_cov);
        ll = ll-log(det(err_cov))/2-err'*inv_err_cov*err/2;
        grad_err = zeros(J,length(param));
        for i=1:length(param)
            grad_err(:,i) = -H*grad_x_pred1(:,i)-grad_H(:,:,i)*x_pred1;
        end
        grad_err_cov = zeros(J,J,length(param));
        for i=1:length(param)
            grad_err_cov(:,:,i) = H*grad_V_pred1(:,:,i)*H'+grad_H(:,:,i)*V_pred1*H'+H*V_pred1*grad_H(:,:,i)'+grad_R(:,:,i);
        end
        hess_err = zeros(J,length(param),length(param));
        hess_err_cov = zeros(J,J,length(param),length(param));
        for i=1:length(param)
            for j=i:length(param)
                hess_err(:,i,j) = -grad_H(:,:,i)*grad_x_pred1(:,j)-grad_H(:,:,j)*grad_x_pred1(:,i)-H*hess_x_pred1(:,i,j);
                hess_err(:,j,i) = hess_err(:,i,j);
                hess_err_cov(:,:,i,j) = grad_H(:,:,i)*grad_V_pred1(:,:,j)*H'+grad_H(:,:,j)*grad_V_pred1(:,:,i)*H'+H*hess_V_pred1(:,:,i,j)*H'+H*grad_V_pred1(:,:,i)*grad_H(:,:,j)'+H*grad_V_pred1(:,:,j)*grad_H(:,:,i)'+grad_H(:,:,i)*V_pred1*grad_H(:,:,j)'+grad_H(:,:,j)*V_pred1*grad_H(:,:,i)';
                hess_err_cov(:,:,j,i) = hess_err_cov(:,:,i,j);
            end
        end
        for i=1:length(param)
            grad(i) = grad(i)-(trace(inv_err_cov*grad_err_cov(:,:,i))+2*err'*inv_err_cov*grad_err(:,i)-err'*inv_err_cov*grad_err_cov(:,:,i)*inv_err_cov*err)/2;
        end
        for i=1:length(param)
            for j=i:length(param)
                hess(i,j) = hess(i,j)-(-trace(inv_err_cov*grad_err_cov(:,:,i)*inv_err_cov*grad_err_cov(:,:,j))+trace(inv_err_cov*hess_err_cov(:,:,i,j)))/2;
                hess(i,j) = hess(i,j)-(2*grad_err(:,i)'*inv_err_cov*grad_err(:,j)-2*err'*inv_err_cov*grad_err_cov(:,:,j)*inv_err_cov*grad_err(:,i)+2*err'*inv_err_cov*hess_err(:,i,j))/2;
                hess(i,j) = hess(i,j)-(-2*grad_err(:,j)'*inv_err_cov*grad_err_cov(:,:,i)*inv_err_cov*err+2*err'*inv_err_cov*grad_err_cov(:,:,j)*inv_err_cov*grad_err_cov(:,:,i)*inv_err_cov*err-err'*inv_err_cov*hess_err_cov(:,:,i,j)*inv_err_cov*err)/2;
                hess(j,i) = hess(i,j);
            end
        end
        if t == T
            break
        end
        Kg = V_pred1*H'*inv_err_cov;
        grad_Kg = zeros(2*K,J,length(param));
        for i=1:length(param)
            grad_Kg(:,:,i) = (grad_V_pred1(:,:,i)*H'+V_pred1*grad_H(:,:,i)')*inv_err_cov - V_pred1*H'*inv_err_cov*grad_err_cov(:,:,i)*inv_err_cov;
        end
        hess_Kg = zeros(2*K,J,length(param),length(param));
        for i=1:length(param)
            for j=i:length(param)
                hess_Kg(:,:,i,j) = (hess_V_pred1(:,:,i,j)*H'+grad_V_pred1(:,:,i)*grad_H(:,:,j)'+grad_V_pred1(:,:,j)*grad_H(:,:,i)')*inv_err_cov;
                hess_Kg(:,:,i,j) = hess_Kg(:,:,i,j)-(grad_V_pred1(:,:,i)*H'+V_pred1*grad_H(:,:,i)')*inv_err_cov*grad_err_cov(:,:,j)*inv_err_cov;
                hess_Kg(:,:,i,j) = hess_Kg(:,:,i,j)-(grad_V_pred1(:,:,j)*H'+V_pred1*grad_H(:,:,j)')*inv_err_cov*grad_err_cov(:,:,i)*inv_err_cov;
                hess_Kg(:,:,i,j) = hess_Kg(:,:,i,j)+V_pred1*H'*(inv_err_cov*grad_err_cov(:,:,j)*inv_err_cov*grad_err_cov(:,:,i)*inv_err_cov-inv_err_cov*hess_err_cov(:,:,i,j)*inv_err_cov+inv_err_cov*grad_err_cov(:,:,i)*inv_err_cov*grad_err_cov(:,:,j)*inv_err_cov);
                hess_Kg(:,:,j,i) = hess_Kg(:,:,i,j);
            end
        end
        x_filt = x_pred1 + Kg*err;
        V_filt = (eye(2*K)-Kg*H)*V_pred1;
        for i=1:length(param)
            grad_x_filt(:,i) = grad_x_pred1(:,i) + Kg*grad_err(:,i) + grad_Kg(:,:,i)*err;
            grad_V_filt(:,:,i) = grad_V_pred1(:,:,i) - grad_Kg(:,:,i)*H*V_pred1 - Kg*grad_H(:,:,i)*V_pred1 - Kg*H*grad_V_pred1(:,:,i);
        end
        for i=1:length(param)
            for j=i:length(param)
                hess_x_filt(:,i,j) = hess_x_pred1(:,i,j) + grad_Kg(:,:,i)*grad_err(:,j) + grad_Kg(:,:,j)*grad_err(:,i) + Kg*hess_err(:,i,j) + hess_Kg(:,:,i,j)*err;
                hess_x_filt(:,j,i) = hess_x_filt(:,i,j);
                hess_V_filt(:,:,i,j) = hess_V_pred1(:,:,i,j) - hess_Kg(:,:,i,j)*H*V_pred1 - grad_Kg(:,:,i)*grad_H(:,:,j)*V_pred1 - grad_Kg(:,:,i)*H*grad_V_pred1(:,:,j) - grad_Kg(:,:,j)*grad_H(:,:,i)*V_pred1 - Kg*grad_H(:,:,i)*grad_V_pred1(:,:,j) - Kg*grad_H(:,:,i)*grad_V_pred1(:,:,j) - grad_Kg(:,:,j)*H*grad_V_pred1(:,:,i) - Kg*grad_H(:,:,j)*grad_V_pred1(:,:,i) - Kg*H*hess_V_pred1(:,:,i,j);
                hess_V_filt(:,:,j,i) = hess_V_filt(:,:,i,j);
            end
        end
        x_pred1 = F*x_filt;
        V_pred1 = F*V_filt*F'+Q;
        for i=1:length(param)
            grad_x_pred1(:,i) = F*grad_x_filt(:,i)+grad_F(:,:,i)*x_filt;
            grad_V_pred1(:,:,i) = F*grad_V_filt(:,:,i)*F' + grad_F(:,:,i)*V_filt*F' + F*V_filt*grad_F(:,:,i)' + grad_Q(:,:,i);
        end
        for i=1:length(param)
            for j=i:length(param)
                hess_x_pred1(:,i,j) = grad_F(:,:,i)*grad_x_filt(:,j) + grad_F(:,:,j)*grad_x_filt(:,i) + F*hess_x_filt(:,i,j) + hess_F(:,:,i,j)*x_filt;
                hess_x_pred1(:,j,i) = hess_x_pred1(:,i,j);
                hess_V_pred1(:,:,i,j) = grad_F(:,:,i)*grad_V_filt(:,:,j)*F' + grad_F(:,:,j)*grad_V_filt(:,:,i)*F' + F*hess_V_filt(:,:,i,j)*F' + F*grad_V_filt(:,:,i)*grad_F(:,:,j)' + F*grad_V_filt(:,:,j)*grad_F(:,:,i)' + hess_F(:,:,i,j)*V_pred1*F' + grad_F(:,:,i)*V_pred1*grad_F(:,:,j)'+ grad_F(:,:,j)*V_pred1*grad_F(:,:,i)' + F*V_pred1*hess_F(:,:,i,j)';
                hess_V_pred1(:,:,j,i) = hess_V_pred1(:,:,i,j);
            end
        end
    end
    grad(K+1:2*K) = grad(K+1:2*K)*2*pi/fs;
    hess(K+1:2*K,:) = hess(K+1:2*K,:)*2*pi/fs;
    hess(:,K+1:2*K) = hess(:,K+1:2*K)*2*pi/fs;
    
    mll = -ll;
    grad = -grad;
    hess = -hess;
end


