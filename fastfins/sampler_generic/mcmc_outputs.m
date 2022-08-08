function out = mcmc_outputs(np, options)
%Create MCMC outputs
%
%Tiangang Cui, August, 2019

m = floor(options.nstep/options.sbatch);

out.mlpt    = zeros(1,m);
out.mllkd   = zeros(1,m);
out.sigma   = zeros(1,m);
out.mh      = zeros(1,m);
out.samples = zeros(np,m);
out.j       = 0;

end