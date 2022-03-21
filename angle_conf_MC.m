function conf = angle_conf_MC(mu,Sigma,prob,seeds)
    nmc = size(seeds,2);
    phi = atan2(mu(2),mu(1));
    P = [cos(phi) sin(phi); -sin(phi) cos(phi)];
    Sigmasqrt = sqrtm(Sigma);
    seeds = norm(mu)*[ones(1,nmc); zeros(1,nmc)]+P*Sigmasqrt*seeds;
    phases = atan2(seeds(2,:),seeds(1,:));
    [tmp,I] = sort(abs(phases));
    phases = sort(phases(I(1:ceil(prob*nmc))));
    conf = [phi+phases(1) phi+phases(length(phases))];
    if conf(1)<-pi
        conf(1) = conf(1)+2*pi;
    end
    if conf(2)>pi
        conf(2) = conf(2)-2*pi;
    end
end

