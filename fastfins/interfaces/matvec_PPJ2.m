function J2dv = matvec_PPJ2(model, prior, sol, dv)
%MATVEC_PPH
%
% compute the matvec with PPH
%
% Tiangang Cui, 19/Mar/2014

du = matvec_prior_L(prior, dv); % transform to u
J2dv = matvec_J2du(model, sol, du);
for i = 1:length(J2dv)
    J2dv{i} = matvec_prior_Lt(prior, J2dv{i}')';
end

end
