function w = matvec_PPFisher(model, prior, sol, dv)
%MATVEC_PPH
%
% compute the matvec with PPH
%
% Tiangang Cui, 19/Mar/2014

du = matvec_prior_L(prior, dv); % transform to u
wu = zeros(size(du));
% matvec with GNH for the model
for i = 1:size(du,2)
    tmp     = sol.I*matvec_Ju (model, sol, du(:,i));
    wu(:,i) = matvec_Jty(model, sol, tmp);
end
w = matvec_prior_Lt(prior, wu); % transform back to v

end
