function [VARwithnoise_A,VARwithnoise_E,VARwithnoise_r,VARwithnoise_AIC] = VAR_fit(Y,MAX_VAR)
	J = size(Y,1);
    VARwithnoise_AIC = zeros(1,MAX_VAR);
    VARwithnoise_A = zeros(J,J*MAX_VAR,MAX_VAR);
    VARwithnoise_E = zeros(J,J,MAX_VAR);
    VARwithnoise_r = zeros(1,MAX_VAR);
    for ARdeg=1:MAX_VAR
        [A,E,r,mll] = VAR_myule(Y,ARdeg);
        VARwithnoise_A(:,1:J*ARdeg,ARdeg) = -A(:,J+1:(ARdeg+1)*J);
        VARwithnoise_E(:,:,ARdeg) = E;
        VARwithnoise_r(ARdeg) = r;
        VARwithnoise_AIC(ARdeg) = 2*mll+2*(J^2*ARdeg+J*(J+1)/2+1);
    end
end
