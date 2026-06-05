function [A,E,R] = armyule(y,ARdeg)
    y = y(:)';
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
    
    R = fminbnd(@(R_val) obj_fun(R_val, y, C, c, acov, ARdeg), 0, min(eigs), options);
    
    A = (C-R*eye(ARdeg))\c;
    E = acov(1)-R-A'*acov(2:ARdeg+1);
    A = [1 -A'];
end

function mll = obj_fun(R_val, y, C, c, acov, ARdeg)
    A_temp = (C - R_val * eye(ARdeg)) \ c;
    log_E = log(acov(1) - R_val - A_temp' * acov(2:ARdeg+1));
    param = [-A_temp; log_E - log(R_val)]';
    mll = ARwithnoise_ll(y, param);
end