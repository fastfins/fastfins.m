function [mllkd_pm, mllkd, v, gmllkd_pm, gmllkd]= minus_sub_pm(mlpost, sub_V, beta, v_sub, pm_size)
%minus_sub_PM
%
% computes the pseudo marginal
%
% Tiangang Cui, 19/Mar/2020

[m,n] = size(v_sub);
mllkd_pm = zeros(1, n);
%
v = zeros(size(sub_V,1), n);
mllkd = zeros(1, n);
%
if nargout > 3
    gmllkd_pm = zeros(m, n);
    gmllkd = zeros(size(sub_V,1), n);
end
%
ml = zeros(1, pm_size);
mg = zeros(size(sub_V, 1), pm_size);
%
for i = 1:n
    rs      = randn(size(sub_V,1), pm_size);
    v_null  = rs - sub_V*(sub_V'*rs);
    vs      = sub_V*v_sub(:,i) + v_null;
    %
    if nargout > 3
        for j = 1:pm_size
            [~,ml(j),~,mg(:,j)] = mlpost(vs(:,j));
        end
    else
        for j = 1:pm_size
            [~,ml(j)] = mlpost(vs(:,j));
        end
    end
    %
    ref     = max(-ml*beta);
    tmp     = sum(exp(-ml*beta-ref));
    weights = exp(-ml*beta-ref)/tmp;
    % 
    mllkd_pm(i) = (-log(tmp/pm_size)-ref);
    %
    % sample an index
    ind         = datasample(1:pm_size, 1, 'weights', weights);
    v(:,i)      = vs(:,ind);
    mllkd(i)    = ml(ind);
    %
    if nargout > 3
        gmllkd_pm(:,i)  = ( sub_V'*(mg*weights(:)) )*beta;
        gmllkd(:,i) = mg(:,ind);
    end
end

end
