function [mll,A,E] = VAR_myule_ll(Y,r,C,c)
    d = size(Y,1);
    ARdeg = size(C,1)/d;
    
    r_eye = r * eye(size(C,1));
    B = (C - r_eye) \ c;
    
    A = reshape(B', d, d * ARdeg);
    
    E = C(1:d,1:d) - r*eye(d);
    for k=1:ARdeg
        idx = (k-1)*d+1:k*d;
        E = E - A(:,idx) * c(idx,:);
    end
    
    [~,z] = polyeig_VAR(A);
    if max(abs(z)) < 1 && min(eig(E)) > 0
        mll = VARwithnoise_ll(Y,A,E,r*eye(d));
    else
        mll = inf;
    end
end
