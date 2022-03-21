function [A,E,R] = armyule(y,ARdeg)
    T = length(y);
    acov = zeros(ARdeg+1,1);
    for k=1:ARdeg+1
        acov(k) = y(1:T-k+1)*y(k:T)'/T;
    end
    C = zeros(ARdeg,ARdeg);
    for k=1:ARdeg
        C(k,1:k) = acov(k:-1:1);
        C(k,k+1:ARdeg) = acov(2:ARdeg-k+1);
    end
    c = acov(2:ARdeg+1);
    eigs = eig([C flipud(c); fliplr(c') acov(1)]);
    options = optimset('Display','off');
    R = fminbnd(@(R)ARwithnoise_ll(y,[-(C-R*eye(ARdeg))\c; log(acov(1)-R-((C-R*eye(ARdeg))\c)'*acov(2:ARdeg+1))-log(R)]'),0,min(eigs),options);
    A = (C-R*eye(ARdeg))\c;
    E = acov(1)-R-A'*acov(2:ARdeg+1);
    A = [1 -A'];
end

