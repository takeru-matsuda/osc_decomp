function [A,E,r,mll] = VAR_myule(Y,ARdeg)
    d = size(Y,1);
    T = size(Y,2);
    acov = zeros(d,d,ARdeg+1);
    Y = Y-mean(Y,2)*ones(1,T);
    for k=1:ARdeg+1
        acov(:,:,k) = Y(:,1:T-k+1)*Y(:,k:T)'/T;
    end
    C = zeros(d*ARdeg,d*ARdeg);
    for k=1:ARdeg
        for j=1:k
            C((k-1)*d+1:k*d,(j-1)*d+1:j*d) = acov(:,:,k-j+1);
        end
        for j=k+1:ARdeg
            C((k-1)*d+1:k*d,(j-1)*d+1:j*d) = acov(:,:,j-k+1)';
        end
    end
    c = zeros(d*ARdeg,d);
    for k=1:ARdeg
        c((k-1)*d+1:k*d,:) = acov(:,:,k+1);
    end
    [r,mll] = fminbnd(@(r)VAR_myule_ll(Y,r,C,c),0,min(eig(C)));
%    r = 0;
    [mll,A,E] = VAR_myule_ll(Y,r,C,c);
    A = [eye(d) -A];
end

