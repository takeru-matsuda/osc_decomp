function beta = whittle_multi_fit(X,y)
    options = optimoptions('fminunc','Algorithm','quasi-newton','Display','off');
    beta = exp(fminunc(@(b)whittle_multi_ll(X,y,b),zeros(size(X,4),1),options));
end

