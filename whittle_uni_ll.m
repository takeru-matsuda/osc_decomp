function mll = whittle_uni_ll(X,y,beta)
    mll = sum(log(X*beta)+y./(X*beta));
%    grad = zeros(size(X,2),1);
%    for i=1:length(y)
%        grad = grad+X(i,:)'/(X(i,:)*beta)-y(i)*X(i,:)'/(X(i,:)*beta)^2;
%    end
end

