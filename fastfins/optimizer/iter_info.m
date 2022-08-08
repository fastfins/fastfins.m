function info = iter_info(opts, f0, jump, fdiff, fnew, n_gnew, i, iter)
% ITER_INFO  determines the terminition condition
%
% Tiangang Cui, 16/August/2020

info = 0;

% display iteration info
fprintf('%4i \t\t %10.5E \t\t %10.5E \t\t %4i\n', [i, fnew, n_gnew, iter]);

if n_gnew < opts.first_KKT_tol
    info = 1;
    disp('First KKT condition satisfied, exit');
    return;
end

if jump < opts.jump_size_tol
    info = 2;
    disp('Jump size reach minumum threshold, exit')
    return
end

if fdiff/abs(f0) < opts.fval_tol
    info = 3;
    disp('Objective function does not have suffcient decrease, exit');
    return
end

if i > opts.max_iter
    info = 4;
    disp('Maximum number of iteration reached, no optimal solution find');
    return;
end

end