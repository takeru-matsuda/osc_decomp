% slight modification of polyeig
function [X,E,s] = mypolyeig(C)
    n = length(C{1});
    p = length(C)-1;
    if nargout > 2 && p < 1
        error('MATLAB:polyeig:tooFewInputs', ...
          'Must provide at least two matrices.')
    end

    A = eye(n*p);
    A(1:n,1:n) = C{1};
    if p == 0
        B = eye(n);
        p = 1;
    else
        B = diag(ones(n*(p-1),1),-n);
        j = 1:n;
        for k = 1:p
            B(1:n,j) = - C{k+1};
            j = j+n;
        end
    end

    % Use the QZ algorithm on the big matrix pair (A,B).
    if nargout > 1
        [X,E] = eig(A,B);
        E = diag(E);
    else
        X = eig(A,B);
        return
    end

    if p >= 2

        % For each eigenvalue, extract the eigenvector from whichever portion
        % of the big eigenvector matrix X gives the smallest normalized residual.
        V = zeros(n,p);
        for j = 1:p*n
            V(:) = X(:,j);
            R = C{p+1};
            if ~isinf(E(j))
                for k = p:-1:1
                    R = C{k} + E(j)*R;
                end
            end
            R = R*V;
           res = sum(abs(R))./ sum(abs(V));  % Normalized residuals.
            [~,ind] = min(res);
            X(1:n,j) = V(:,ind)/norm(V(:,ind));  % Eigenvector with unit 2-norm.
        end
        X = X(1:n,:);

    end

    if nargout > 2
        % Construct matrix Y whose rows are conjugates of left eigenvectors.
        rcond_p = rcond(C{p+1});
        rcond_0 = rcond(C{1});
        if max(rcond_p,rcond_0) <= eps
            error('MATLAB:polyeig:nonSingularCoeffMatrix', ...
              'Either the leading or the trailing coefficient matrix must be nonsingular.')
        end
        if rcond_p >= rcond_0
            V = C{p+1}; E1 = E;
        else
            V = C{1}; E1 = 1./E;
        end
        Y = X;
        for i=1:p-1
            Y = [Y; Y(end-n+1:end,:)*diag(E1)];
        end
        B = zeros(p*n,n); B(end-n+1:end,1:n) = eye(n);
        Y = Y\B/V;
        for i = 1:n*p, Y(i,:) = Y(i,:)/norm(Y(i,:)); end % Normalize left eigenvectors.
    
        % Reconstruct alpha-beta representation of eigenvalues: E(i) = a(i)/b(i).
        a = E;
        b = ones(size(a));
        k = isinf(a); a(k) = 1; b(k) = 0;
    
        nu = zeros(p+1,1);
        for i = 1:p+1, nu(i) = norm(C{i},'fro'); end
        s = zeros(n*p,1);
    
        % Compute condition numbers.
        for j = 1:p*n
            ab = (a(j).^(0:p-1)).*(b(j).^(p-1:-1:0));
            Da = ab(1)*C{2};
            Db = p*ab(1)*C{1};
            for k = 2:p
                Da = Da + k*ab(k)*C{k+1};
                Db = Db + (p-k+1)*ab(k)*C{k};
            end
            nab = norm( (a(j).^(0:p)) .* (b(j).^(p:-1:0)) .* nu' );
            s(j) = nab / abs( Y(j,:) * (conj(b(j))*Da-conj(a(j))*Db) * X(:,j) );
        end
    
    end
end
