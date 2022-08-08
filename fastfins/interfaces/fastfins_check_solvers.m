function fastfins_check_solvers()

if exist('forward_solve') == 0
    error('The function forward_solve() does not exist, cannot do anything')
end

if exist('forward_solve_vec') == 0
    disp('The function forward_solve_vec() does not exist, cannot run minus_log_post_vec')
end

if exist('prior_cov_l') == 0
    error('The function prior_cov_l() does not exist, cannot perform forward parameter transformation for forward model simulation')
end

if exist('matvec_Jty') == 0
    warning('The function matvec_Jty() does not exist, cannot compute gradient')
end

if exist('prior_cov_lt') == 0
    warning('The function prior_cov_lt() does not exist, cannot perform back transformation of gradient')
end

if exist('matvec_Ju') == 0
    warning('The function matvec_Ju() does not exist, cannot compute Fisher information')
end

if exist('prior_cov_ilt') == 0
    warning('The function prior_cov_ilt() does not exist, may impact LIS dimension reduction')
end

if exist('prior_cov_c') == 0
    warning('The function prior_cov_c() does not exist, cannot carry KL dimension reduction')
end

if exist('prior_cov_eig') == 0
    warning('The function prior_cov_eig() does not exist, cannot carry KL dimension reduction')
end

if exist('prior_cov_il') == 0
    warning('The function prior_cov_il() does not exist, may impact interpretation of samples and KL dimension reduction')
end

end