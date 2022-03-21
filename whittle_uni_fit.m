function beta = whittle_uni_fit(X,y)
    options = optimoptions('fminunc','Algorithm','quasi-newton','MaxFunEvals',10000,'Display','off');
    beta = exp(fminunc(@(b)whittle_uni_ll(X,y,exp(b)),zeros(size(X,2),1),options));
end

