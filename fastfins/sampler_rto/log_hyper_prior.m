function log_p0 = log_hyper_prior(hyper, thetas)

[m, n] = size(thetas);
log_p0 = zeros(1, n);

for i = 1:n
    for j = 1:m
        switch hyper{j}.type
            case{'gamma'}
                log_p0(i) = log_p0(i) + (hyper{j}.alpha-1)*log(thetas(j,i)) - hyper{j}.beta*thetas(j,i);
            case{'beta'}
                log_p0(i) = log_p0(i) +  hyper{j}.alpha*log(abs(thetas(j,i) - hyper{j}.left)) ...
                    + hyper{j}.beta*log(abs(hyper{j}.right - thetas(j,i)));
        end
    end
end

end

