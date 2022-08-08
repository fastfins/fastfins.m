function [next_h, lr, inc] = hyper_propose(obs, prior, hyper_setup, adapt, curr_h)

c_tmp   = log([curr_h.alpha; curr_h.beta]);

if rand < adapt.AM_weight
    n_tmp   = c_tmp + adapt.L*randn(size(c_tmp));
    lr      = 0;
    inc     = 1;
    disp('AM')
else
    n_tmp   = random(adapt.GM, size(c_tmp, 2))';
    %lr      = sum(c_tmp, 1)' - sum(n_tmp, 1)' + log(pdf(adapt.GM, c_tmp')) - log(pdf(adapt.GM, n_tmp'));
    lr      = log(pdf(adapt.GM, c_tmp')) - log(pdf(adapt.GM, n_tmp'));
    inc     = 0;
    disp('GM')
end

n_alpha = exp(n_tmp(hyper_setup.alpha_ind));
n_beta  = exp(n_tmp(hyper_setup.beta_ind));
next_h  = update_model(obs, prior, hyper_setup, n_alpha, n_beta);

end
