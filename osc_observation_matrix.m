function H = osc_observation_matrix(c,J,K)
% c stores adjacent coordinate pairs, first by channel, then by oscillator.
H=zeros(J,2*K);
H(1,1:2:end)=1;
if J==1, return; end
assert(numel(c)==2*(J-1)*K,'osc_decomp:LoadingSize', ...
    'Expected two observation coefficients per channel and oscillator.');
for k=1:K
    block=c((k-1)*2*(J-1)+(1:2*(J-1)));
    H(2:end,2*k-1:2*k)=reshape(block,2,J-1)';
end
end
