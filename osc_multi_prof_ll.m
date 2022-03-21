% profile likelihood (observation noise variance reduced by Kitagawa method)
function [mll,grad,tau2hat] = osc_multi_prof_ll(Y,param,init_theta,return_grad)
    J = size(Y,1);
    T = size(Y,2);
    K = length(param)/(2*J+1);
    a = (tanh(param(1:K))+1)/2;
    theta = init_theta+tanh(param(K+1:2*K))*pi;
    sigma2 = exp(param(2*K+1:3*K));
    c = param(3*K+1:end);
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
    R = eye(J);
    
    x_pred1 = zeros(2*K,T);
    x_filt = zeros(2*K,T);
    V_pred1 = zeros(2*K,2*K,T);
    V_filt = zeros(2*K,2*K,T);
    x_pred1(:,1) = zeros(2*K,1);
    for k=1:K
        V_pred1(2*k-1:2*k,2*k-1:2*k,1) = sigma2(k)/(1-a(k)^2)*eye(2);
    end
    for t=1:T-1
        x_filt(:,t) = x_pred1(:,t) + V_pred1(:,:,t)*H'*((H*V_pred1(:,:,t)*H'+R)\(Y(:,t)-H*x_pred1(:,t)));
        V_filt(:,:,t) = V_pred1(:,:,t) - V_pred1(:,:,t)*H'*((H*V_pred1(:,:,t)*H'+R)\H)*V_pred1(:,:,t);
        x_pred1(:,t+1) = F*x_filt(:,t);
        V_pred1(:,:,t+1) = F*V_filt(:,:,t)*F'+Q;
    end
    tau2hat = 0;
    for t=1:T
        tau2hat = tau2hat+(Y(:,t)-H*x_pred1(:,t))'*((H*V_pred1(:,:,t)*H'+R)\(Y(:,t)-H*x_pred1(:,t)))/J/T;
    end
    ll = -J*T*log(tau2hat)/2-J*T/2-J*T/2*log(2*pi);
    for t=1:T
        ll = ll-log(det(H*V_pred1(:,:,t)*H'+R))/2;
    end
    mll = -ll;
    
    if return_grad == false
        grad = 0;
        return
    end
    grad_F = zeros(2*K,2*K,(2*J+1)*K);
    grad_Q = zeros(2*K,2*K,(2*J+1)*K);
    grad_H = zeros(J,2*K,(2*J+1)*K);
    kk = 1;
    for k=1:K
        grad_F(2*k-1:2*k,2*k-1:2*k,k) = 1/2/cosh(param(k))^2*[cos(theta(k)) -sin(theta(k)); sin(theta(k)) cos(theta(k))];
        grad_F(2*k-1:2*k,2*k-1:2*k,K+k) = pi/cosh(param(K+k))^2*a(k)*[-sin(theta(k)) -cos(theta(k)); cos(theta(k)) -sin(theta(k))];
        grad_Q(2*k-1:2*k,2*k-1:2*k,2*K+k) = sigma2(k)*eye(2);
        for j=2:J
            grad_H(j,2*k-1,3*K+kk) = 1;
            grad_H(j,2*k,3*K+kk+1) = 1;
            kk = kk+2;
        end
    end
    grad_x_pred1 = zeros(2*K,(2*J+1)*K);
    grad_V_pred1 = zeros(2*K,2*K,(2*J+1)*K);
    for k=1:K
        grad_V_pred1(2*k-1:2*k,2*k-1:2*k,k) = a(k)/cosh(param(k))^2*sigma2(k)/(1-a(k)^2)^2*eye(2);
        grad_V_pred1(2*k-1:2*k,2*k-1:2*k,2*K+k) = sigma2(k)/(1-a(k)^2)*eye(2);
    end
    grad_x_filt = zeros(2*K,(2*J+1)*K);
    grad_V_filt = zeros(2*K,2*K,(2*J+1)*K);
    grad_tau2hat = zeros((2*J+1)*K,1);
    grad = zeros((2*J+1)*K,1);
    for t=1:T
        err = Y(:,t)-H*x_pred1(:,t);
        err_cov = H*V_pred1(:,:,t)*H'+R;
        inv_err_cov = inv(err_cov);
        grad_err = zeros(J,(2*J+1)*K);
        for i=1:(2*J+1)*K
            grad_err(:,i) = -H*grad_x_pred1(:,i)-grad_H(:,:,i)*x_pred1(:,t);
        end
        grad_err_cov = zeros(J,J,(2*J+1)*K);
        for i=1:(2*J+1)*K
            grad_err_cov(:,:,i) = H*grad_V_pred1(:,:,i)*H'+grad_H(:,:,i)*V_pred1(:,:,t)*H'+H*V_pred1(:,:,t)*grad_H(:,:,i)';
        end
        for i=1:(2*J+1)*K
            grad(i) = grad(i)-(trace(inv_err_cov*grad_err_cov(:,:,i)))/2;
            grad_tau2hat(i) = grad_tau2hat(i)-2*(Y(:,t)-H*x_pred1(:,t))'*inv_err_cov*(H*grad_x_pred1(:,i)+grad_H(:,:,i)*x_pred1(:,t))/J/T;
            grad_tau2hat(i) = grad_tau2hat(i)-(Y(:,t)-H*x_pred1(:,t))'*inv_err_cov*grad_err_cov(:,:,i)*inv_err_cov*(Y(:,t)-H*x_pred1(:,t))/J/T;
        end
        if t == T
            break
        end
        Kg = V_pred1(:,:,t)*H'*inv_err_cov;
        grad_Kg = zeros(2*K,J,(2*J+1)*K);
        for i=1:(2*J+1)*K
            grad_Kg(:,:,i) = (grad_V_pred1(:,:,i)*H'+V_pred1(:,:,t)*grad_H(:,:,i)')*inv_err_cov - V_pred1(:,:,t)*H'*inv_err_cov*grad_err_cov(:,:,i)*inv_err_cov;
        end
        for i=1:(2*J+1)*K
            grad_x_filt(:,i) = grad_x_pred1(:,i) + Kg*grad_err(:,i) + grad_Kg(:,:,i)*err;
            grad_V_filt(:,:,i) = grad_V_pred1(:,:,i) - grad_Kg(:,:,i)*H*V_pred1(:,:,t) - Kg*grad_H(:,:,i)*V_pred1(:,:,t) - Kg*H*grad_V_pred1(:,:,i);
        end
        for i=1:(2*J+1)*K
            grad_x_pred1(:,i) = F*grad_x_filt(:,i)+grad_F(:,:,i)*x_filt(:,t);
            grad_V_pred1(:,:,i) = F*grad_V_filt(:,:,i)*F' + grad_F(:,:,i)*V_filt(:,:,t)*F' + F*V_filt(:,:,t)*grad_F(:,:,i)' + grad_Q(:,:,i);
        end
    end
    grad = grad-J*T*grad_tau2hat/tau2hat/2;
    grad = -grad;
end


