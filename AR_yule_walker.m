function [A,E] = AR_yule_walker(y,p)
    arguments
        y {mustBeNumeric,mustBeVector,mustBeNonempty,mustBeFinite}
        p (1,1) {mustBeInteger,mustBeNonnegative}
    end
    y = double(y(:));
    T = numel(y);
    if p >= T
        error('AR_yule_walker:OrderTooLarge','The AR order must be smaller than the sample count.');
    end
    r = zeros(p+1,1);
    for lag = 0:p
        r(lag+1) = y(1:T-lag)'*y(1+lag:T)/T;
    end
    A = 1;
    E = real(r(1));
    if ~(isfinite(E) && E > 0)
        error('AR_yule_walker:InvalidVariance','The signal must have finite, positive mean-square amplitude.');
    end
    for order = 1:p
        reflection = -(r(order+1)+A(2:order)*r(order:-1:2))/E;
        A = [A 0]+reflection*[0 conj(fliplr(A))];
        E = E*(1-abs(reflection)^2);
        if ~(isfinite(E) && E > 0)
            error('AR_yule_walker:InvalidVariance','The estimated prediction variance must be positive.');
        end
    end
end
