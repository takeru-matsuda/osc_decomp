function [A,E] = AR_MLE(y,ARdeg)
    A = aryule(y,ARdeg);
    a = zeros(ARdeg,ARdeg);
    a(:,ARdeg) = -A(2:ARdeg+1);
    for m=ARdeg:-1:2
        for i=1:m-1
            a(i,m-1) = (a(i,m)+a(m,m)*a(m-i,m))/(1-a(m,m)^2);
        end
    end
    c = zeros(ARdeg,1); % PARCOR
    for m=1:ARdeg
        c(m) = a(m,m);
    end
    init = log(1+c)-log(1-c);
    init(init>20) = 20;
    init(init<-20) = -20;
    init(c>=1) = 20;
    init(c<=-1) = -20;
    options = optimoptions('fminunc','Display','off');
    tmp = fminunc(@(p)AR_ll(y,p),init,options);
    c = (exp(tmp)-1)./(exp(tmp)+1);
    [~,E] = AR_ll(y,tmp);
    for m=1:ARdeg
        a(m,m) = c(m);
    end
    for m=2:ARdeg
        for i=1:m-1
            a(i,m) = a(i,m-1)-c(m)*a(m-i,m-1);
        end
    end
    A = [1 -a(:,ARdeg)'];
end

