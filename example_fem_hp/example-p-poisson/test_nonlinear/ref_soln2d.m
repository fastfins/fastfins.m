function [state_func, eta_func, force_func, beta_func] = ref_soln2d(p_rate, epsilon, plot_on)
%ref_soln2d
% 
% generates reference solution and corresponding forcing 
%
% Tiangang Cui, 17/Oct/2016

syms x1 x2

%state_func(x1, x2) = ( 1 - cos(2*pi*x1) ) * sin( 1.25*pi*x2 + 0.25*pi );
state_func(x1, x2) = ( 1 - cos(2*pi*x1) ) * cos(0.25*pi*x2 - 0.25*pi);

ds_x1(x1, x2) = diff(state_func, x1);
ds_x2(x1, x2) = diff(state_func, x2);

% eta
grad_s_l2n(x1, x2) = ds_x1(x1, x2)^2 + ds_x2(x1, x2)^2;
eta_func(x1, x2) = (grad_s_l2n(x1, x2) + epsilon)^(p_rate/2-1);

% forcing
force_func(x1, x2) = - (diff(eta_func*ds_x1, x1) + ...
    diff(eta_func*ds_x2, x2));

% robin at bottom boundary
grad_s_dot_n(x1) = -ds_x2(x1, 0);
beta_func(x1) = -eta_func(x1,0)*grad_s_dot_n(x1)/state_func(x1,0);

% eval(int(state_func*force_func, x1, 0, 1))


% check left and right boundaries
disp('left flux:')
disp(eta_func(0, x2)*ds_x1(0, x2))
disp('right flux:')
disp(eta_func(1, x2)*ds_x1(1, x2))
% check top boundary
disp('top flux:')
disp(eta_func(x1, 1)*ds_x2(x1, 1))

disp('top Dirichlet:')
disp(state_func(x1, 1))

if plot_on
    xs = linspace(0,1,50);
    ys = xs;
    
    figure
    subplot(2,2,1)
    plot(xs, eval(eta_func(xs, 0)))
    title('\eta(x_1, 0)')
    
    subplot(2,2,2)
    plot(xs, eval(grad_s_dot_n(xs)))
    title('grad(x_1, x_2) \cdot n')
    
    subplot(2,2,3)
    plot(xs, eval(state_func(xs,0)))
    title('u(x_1, 0)')
    
    subplot(2,2,4)
    plot(xs, eval(beta_func(xs)))
    title('\beta(x_1)')
    
    [xx, yy] = meshgrid(xs, ys);
    ss = eval(state_func(xx, yy));
    ff = eval(force_func(xx, yy));
    
    figure
    surf(xx, yy, ss)
    xlabel('x_1')
    ylabel('x_2')
    title('u_{ref}(x_1, x_2)')
    
    figure
    surf(xx, yy, ff)
    xlabel('x_1')
    ylabel('x_2')
    title('f(x_1, x_2)')
end


state_func  = matlabFunction(state_func);
force_func  = matlabFunction(force_func);
eta_func    = matlabFunction(eta_func);
beta_func   = matlabFunction(beta_func);

end