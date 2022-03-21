function mll = whittle_multi_ll(X,y,logbeta)
    d = size(X,1);
    K = size(X,4);
    beta = exp(logbeta);
    mll = 0;
%    grad = zeros(K,1);
    for i=1:size(y,2)
        F = zeros(d,d);
        for k=1:K
            F = F+2*pi*X(:,:,i,k)*beta(k);
        end
        mll = mll+log(det(F))+y(:,i)'*(F\y(:,i));
%        for k=1:K
%            grad(k) = grad(k)+(trace(F\X(:,:,i,k))-(F\y(:,i))'*X(:,:,i,k)*(F\y(:,i)))*beta(k);
%        end
    end
    mll = real(mll);
end

