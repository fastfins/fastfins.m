function [state_f, kappa_f, forcing_f, bottom_flux, bnd_funs, bc_types, bc_funs] = ref_soln2d_scalar(plot_on)
%ref_soln2d
% 
% generates reference solution and corresponding forcing 
% bottom and top: Dirichlet
% left and right: no flux
%
% Tiangang Cui, 17/Oct/2016

syms x1 x2

kappa(x1, x2) = 2 + cos(2*pi*x1 + pi/4) * sin(2*pi*x2 + pi/4);

% state_func(x1, x2) = ( 1 - cos(2*pi*x1) ) * sin( 0.25*pi*(1-x2) );
state(x1, x2) = ( 1 - cos(2*pi*x1) ) * sin(1.25*pi*x2 + 0.25*pi );

ds_x1(x1, x2) = diff(state, x1);
ds_x2(x1, x2) = diff(state, x2);

% forcing
forcing(x1, x2) = - (diff(kappa*ds_x1, x1) + diff(kappa*ds_x2, x2));

% flux
left_flux  (x1, x2) = kappa(0, x2)*ds_x1(0, x2);
right_flux (x1, x2) = kappa(1, x2)*ds_x1(1, x2);
bottom_flux(x1, x2) = kappa(x1, 0)*ds_x2(x1, 0);
top_flux   (x1, x2) = kappa(x1, 1)*ds_x2(x1, 1);

left_flux   = matlabFunction(left_flux);
right_flux  = matlabFunction(right_flux);
bottom_flux = matlabFunction(bottom_flux);
top_flux    = matlabFunction(top_flux);

% essential 
left_ess  (x1, x2)  = state(0, x2);
right_ess (x1, x2)  = state(1, x2);
bottom_ess(x1, x2)  = state(x1, 0);
top_ess   (x1, x2)  = state(x1, 1);

left_ess    = matlabFunction(left_ess);
right_ess   = matlabFunction(right_ess);
bottom_ess  = matlabFunction(bottom_ess);
top_ess     = matlabFunction(top_ess);

disp('left flux:')
disp(left_flux)
disp('bottom flux:')
disp(bottom_flux)
disp('right flux:')
disp(right_flux)
disp('top flux:')
disp(top_flux)

disp('left Dirichlet:')
disp(left_ess)
disp('bottom Dirichlet:')
disp(bottom_ess)
disp('right Dirichlet:')
disp(right_ess)
disp('top Dirichlet:')
disp(top_ess)

% left, bottom, right, top
tol = 1E-4;
%bnd_funs = {@(p) find(abs(p(:,1))<tol), @(p) find(abs(p(:,2))<tol), @(p) find(abs(p(:,1)-1)<tol), @(p) find(abs(p(:,2)-1)<tol)};
%bc_types = {'flux', 'essential', 'flux', 'essential'};
%bc_funs  = {left_flux, bottom_ess, right_flux, top_ess};

bnd_funs = {@(p) find(abs(p(:,2))<tol), @(p) find(abs(p(:,2)-1)<tol)};
bc_types = {'essential', 'essential'};
bc_funs  = {bottom_ess, top_ess};

state_f     = matlabFunction(state);
forcing_f   = matlabFunction(forcing);
kappa_f     = matlabFunction(kappa);

if plot_on
    xs = linspace(0,1,50);
    ys = xs;
    
    %{
    figure
    subplot(2,2,1)
    plot(xs, eval(state_fun(xs, 0)))
    title('u(x_1, 0)')
    
    subplot(2,2,2)
    plot(xs, eval(grad_s_dot_n(xs)))
    title('grad(x_1, x_2) \cdot n')
    
    subplot(2,2,3)
    plot(xs, eval(state_func(xs,0)))
    title('u(x_1, 0)')
    
    subplot(2,2,4)
    plot(xs, eval(beta_func(xs)))
    title('\beta(x_1)')
    %}
    
    
    [xx, yy] = meshgrid(xs, ys);
    ss = state_f(xx, yy);
    ff = forcing_f(xx, yy);
    kk = kappa_f(xx, yy);
    figure
    surf(xx, yy, ss)
    xlabel('$x_1$')
    ylabel('$x_2$')
    title('$u_{ref}(x_1, x_2)$')
    
    figure
    surf(xx, yy, ff)
    xlabel('$x_1$')
    ylabel('$x_2$')
    title('$f(x_1, x_2)$')
    
    figure
    surf(xx, yy, kk)
    xlabel('$x_1$')
    ylabel('$x_2$')
    title('$\kappa(x_1, x_2)$')
end

end