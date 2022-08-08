function out = pcn(ml_target, init, options)
%PCN
% Prior preconditioned Crank-Nicolson MCMC
%
% Tiangang Cui, 5/Mar/2013%


curr.v  = init;
np      = length(init);
sigma   = options.sigma;
out     = mcmc_outputs(np, options);

% initialize the MCMC
switch options.proposal
    case {'OW_Prior'}
        [curr.mlpt,curr.mllkd] = ml_target(curr.v);
    case {'MALA'}
        [curr.mlpt,curr.mllkd,~,curr.grad] = ml_target(curr.v);
end

% MCMC
acc_t = 0;
acc   = 0;
batch = 0;

for i = 1:options.nstep
    % propose, and evaluate the acceptance rate
    dt  = exp(sigma);
    a   = 2*sqrt(2*dt)/(2+dt);  % compute beta
    b   = (2-dt)/(2+dt);        % sqrt( 1-beta^2 )
    c   = -2*dt/(2+dt);         % sqrt( 1-beta^2 ) - 1
    
    % propose a condidate, then evaluate the MH ratio
    switch options.proposal
        case {'OW_Prior'}
            r       = randn(np,1);
            next.v  = b*curr.v + a*r;
            % next.lpotential = options.lpotential(next.v);
            [next.mlpt, next.mllkd] = ml_target(next.v);
            alpha   = curr.mllkd - next.mllkd;
        case {'MALA'}
            r       = randn(np,1);
            next.v  = b*curr.v + c*curr.grad + a*r;
            [next.mlpt,next.mllkd,~,next.grad] = ml_target(next.v);
            % reverse
            log_yx  = - next.mllkd - 0.25*dt*norm(next.grad)^2 ...
                - 0.5*next.grad(:)'*(curr.v(:) - next.v(:)) ...
                - 0.25*dt*next.grad(:)'*(curr.v(:) + next.v(:));
            % forward
            log_xy  = - curr.mllkd - 0.25*dt*norm(curr.grad)^2 ...
                - 0.5*curr.grad(:)'*(next.v(:) - curr.v(:)) ...
                - 0.25*dt*curr.grad(:)'*(next.v(:) + curr.v(:));
            alpha   = log_yx - log_xy;
    end
    
    if  log(rand)   < alpha
        curr  = next;
        acc   = acc + 1;
        acc_t = acc_t + 1;
    end
    
    if  options.adapt
        batch = batch + 1;
        if  batch == options.nbatch
            delta = min(0.1,sqrt(options.nbatch/i));
            if (acc/options.nbatch) < options.rate
                sigma = sigma - delta;
            else
                sigma = sigma + delta;
            end
            batch = 0;
            acc   = 0;
        end
    end
    
    if options.sbatch > 1
        if mod(i, options.sbatch) == 0
            out.j = out.j + 1;
            out.samples(:,out.j) = curr.v;
            out.mlpt(out.j)      = curr.mlpt;
            out.mllkd(out.j)     = curr.mllkd;
            out.mh(out.j)        = alpha;
            out.sigma(out.j)     = sigma;
        end
    else
        out.samples(:,i) = curr.v;
        out.mlpt(i)      = curr.mlpt;
        out.mllkd(i)     = curr.mllkd;
        out.mh(i)        = alpha;
        out.sigma(i)     = sigma;
    end
end

out.acc_rate = acc_t/options.nstep;
    
end
