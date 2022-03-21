function [ARwithnoise_param,ARwithnoise_AIC,AR_param,AR_AIC] = AR_fit(y,MAX_AR)
    AR_AIC = zeros(1,MAX_AR);
    AR_param = zeros(MAX_AR,MAX_AR+1);
    ARwithnoise_AIC = zeros(1,MAX_AR);
    ARwithnoise_param = zeros(MAX_AR,MAX_AR+2);
    for ARdeg=1:MAX_AR
%        [A,E] = AR_MLE(y,ARdeg);
        [A,E] = aryule(y,ARdeg);
        a = zeros(ARdeg,ARdeg);
        a(:,ARdeg) = -A(2:ARdeg+1);
        for m=ARdeg:-1:2
            for i=1:m-1
                a(i,m-1) = (a(i,m)+a(m,m)*a(m-i,m))/(1-a(m,m)^2);
            end
        end
        c = zeros(ARdeg,1);
        for m=1:ARdeg
            c(m) = a(m,m);
        end
        AR_AIC(ARdeg) = 2*AR_ll(y,log(1+c)-log(1-c))+2*(ARdeg+1);
        AR_param(ARdeg,1:ARdeg+1) = [A(2:ARdeg+1) E];
        [A,E,R] = armyule(y,ARdeg);
%        param = fminunc(@(param)ARwithnoise_ll(y,param),[A(2:ARdeg+1) log(E)-log(R)],options);
        param = [A(2:ARdeg+1) log(E)-log(R)];
        [mll,R] = ARwithnoise_ll(y,param);
        ARwithnoise_AIC(ARdeg) = 2*mll+2*(ARdeg+2);
        ARwithnoise_param(ARdeg,1:ARdeg+2) = [param(1:ARdeg) exp(param(ARdeg+1))*R R];
    end
end
