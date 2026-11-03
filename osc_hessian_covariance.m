function [covariance,info] = osc_hessian_covariance(hessian)
% Inverse observed information only when it supports a regular Wald interval.
% Return [] with a reason at singular, indefinite, or nonfinite estimates.
covariance=[];
info=struct('valid',false,'reason','','reciprocal_condition',0);
if ~isreal(hessian) || any(~isfinite(hessian),'all')
    info.reason='observed information is not finite and real'; return;
end
H=(hessian+hessian')/2;
[R,flag]=chol(H);
info.reciprocal_condition=rcond(H);
if flag~=0
    info.reason='observed information is not positive definite'; return;
end
if info.reciprocal_condition<1e-12
    info.reason='observed information is numerically singular'; return;
end
covariance=R\(R'\eye(size(H)));
covariance=(covariance+covariance')/2;
info.valid=true;
end
