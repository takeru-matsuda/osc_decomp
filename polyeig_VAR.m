function [V,z] = polyeig_VAR(A)
    d = size(A,1);
    ARdeg = size(A,2)/d;
    C = cell(1,ARdeg+1);
    for k=1:ARdeg
        [C{k}] = -A(:,(ARdeg-k)*d+1:(ARdeg-k+1)*d);
    end
    [C{ARdeg+1}] = eye(d);
    [V,z] = mypolyeig(C);
end

