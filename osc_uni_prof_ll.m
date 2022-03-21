% profile likelihood (observation noise variance reduced by Kitagawa method)
function [mll,grad,tau2hat] = osc_uni_prof_ll(y,param,init_theta,return_grad)
    T = length(y);
    K = length(param)/3;
    a = (tanh(param(1:K))+1)/2;
    theta = init_theta+tanh(param(K+1:2*K))*pi;
    sigma2 = exp(param(2*K+1:3*K));
    F = zeros(2*K,2*K);
    Q = zeros(2*K,2*K);
    for k=1:K
        F(2*k-1:2*k,2*k-1:2*k) = a(k)*[cos(theta(k)) -sin(theta(k)); sin(theta(k)) cos(theta(k))];
        Q(2*k-1:2*k,2*k-1:2*k) = sigma2(k)*eye(2);
    end
    H = zeros(1,2*K);
    H(1:2:2*K) = 1;
    R = 1;

    x_pred1 = zeros(2*K,T);
    x_filt = zeros(2*K,T);
    V_pred1 = zeros(2*K,2*K,T);
    V_filt = zeros(2*K,2*K,T);
    x_pred1(:,1) = zeros(2*K,1);
    for k=1:K
        V_pred1(2*k-1:2*k,2*k-1:2*k,1) = sigma2(k)/(1-a(k)^2)*eye(2);
    end
    for t=1:T-1
        x_filt(:,t) = x_pred1(:,t) + V_pred1(:,:,t)*H'*((H*V_pred1(:,:,t)*H'+R)\(y(t)-H*x_pred1(:,t)));
        V_filt(:,:,t) = V_pred1(:,:,t) - V_pred1(:,:,t)*H'*((H*V_pred1(:,:,t)*H'+R)\H)*V_pred1(:,:,t);
        x_pred1(:,t+1) = F*x_filt(:,t);
        V_pred1(:,:,t+1) = F*V_filt(:,:,t)*F'+Q;
    end
    tau2hat = 0;
    for t=1:T
        tau2hat = tau2hat+(y(t)-H*x_pred1(:,t))^2/(H*V_pred1(:,:,t)*H'+R)/T;
    end
    ll = -T*log(tau2hat)/2-T/2-T/2*log(2*pi);
    for t=1:T
        ll = ll-log(H*V_pred1(:,:,t)*H'+R)/2;
    end
    mll = -ll;

    if return_grad == false
        grad = 0;
        return
    end
    grad_F = zeros(2*K,2*K,3*K);
    grad_Q = zeros(2*K,2*K,3*K);
    for k=1:K
        grad_F(2*k-1:2*k,2*k-1:2*k,k) = 1/2/cosh(param(k))^2*[cos(theta(k)) -sin(theta(k)); sin(theta(k)) cos(theta(k))];
        grad_F(2*k-1:2*k,2*k-1:2*k,K+k) = pi/cosh(param(K+k))^2*a(k)*[-sin(theta(k)) -cos(theta(k)); cos(theta(k)) -sin(theta(k))];
        grad_Q(2*k-1:2*k,2*k-1:2*k,2*K+k) = sigma2(k)*eye(2);
    end
    grad_x_pred1 = zeros(2*K,3*K);
    grad_V_pred1 = zeros(2*K,2*K,3*K);
    for k=1:K
        grad_V_pred1(2*k-1:2*k,2*k-1:2*k,k) = a(k)/cosh(param(k))^2*sigma2(k)/(1-a(k)^2)^2*eye(2);
        grad_V_pred1(2*k-1:2*k,2*k-1:2*k,2*K+k) = sigma2(k)/(1-a(k)^2)*eye(2);
    end
    grad_x_filt = zeros(2*K,3*K);
    grad_V_filt = zeros(2*K,2*K,3*K);
    grad_Rhat = zeros(3*K,1);
    grad = zeros(3*K,1);
    for t=1:T
        err = y(t)-H*x_pred1(:,t);
        err_var = H*V_pred1(:,:,t)*H'+R;
        grad_err = (-H*grad_x_pred1)';
        grad_err_var = zeros(3*K,1);
        for i=1:3*K
            grad_err_var(i) = H*grad_V_pred1(:,:,i)*H';
        end
        for i=1:3*K
            grad(i) = grad(i)-H*grad_V_pred1(:,:,i)*H'/err_var/2;
            grad_Rhat(i) = grad_Rhat(i)-2*(y(t)-H*x_pred1(:,t))*H*grad_x_pred1(:,i)/err_var/T;
            grad_Rhat(i) = grad_Rhat(i)-(y(t)-H*x_pred1(:,t))^2/err_var^2*H*grad_V_pred1(:,:,i)*H'/T;
        end
        if t == T
            break
        end
        Kg = V_pred1(:,:,t)*H'/err_var;
        grad_Kg = zeros(2*K,3*K+1);
        for i=1:3*K
            grad_Kg(:,i) = grad_V_pred1(:,:,i)*H'/err_var-V_pred1(:,:,t)*H'*grad_err_var(i)/err_var^2;
        end
        for i=1:3*K
            grad_x_filt(:,i) = grad_x_pred1(:,i) + Kg*grad_err(i) + grad_Kg(:,i)*err;
            grad_V_filt(:,:,i) = grad_V_pred1(:,:,i) - grad_Kg(:,i)*H*V_pred1(:,:,t) - Kg*H*grad_V_pred1(:,:,i);
        end
        for i=1:3*K
            grad_x_pred1(:,i) = F*grad_x_filt(:,i)+grad_F(:,:,i)*x_filt(:,t);
            grad_V_pred1(:,:,i) = F*grad_V_filt(:,:,i)*F' + grad_F(:,:,i)*V_filt(:,:,t)*F' + F*V_filt(:,:,t)*grad_F(:,:,i)' + grad_Q(:,:,i);
        end
    end
    grad = grad-T*grad_Rhat/tau2hat/2;
    grad = -grad;
end


