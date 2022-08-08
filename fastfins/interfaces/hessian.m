function H = hessian(model, obs, prior, sol)

if obs.jacobian 
    Ju = explicit_jacobian(model, sol);
    Jv = matvec_prior_Lt(prior, Ju')';
    H  = Jv'*sol.I*Jv + eye(prior.dof);
end

end
