clear variables
close all
clc

%% INITIALIZATION

useintervals = false; % change to true if you have Intlab
load('data_SH.mat','u','rho') % A precomputed approximate solution

% Parameters for the proof
if useintervals
    nu = intval('1');
    rstar = intval('1e-5');
else
    nu = 1;
    rstar = 1e-5;
end

% Plot
x = linspace(0, pi, 1e3)';
figure
plot(x, eval_cos(u,x), 'Linewidth', 2)
legend('$\bar{u}(x)$', 'Interpreter', 'Latex', 'Location', 'NorthWest')
xlabel('$x$', 'Interpreter', 'Latex')
title(['Approximate solution for $\rho$ = ',num2str(rho)], 'Interpreter', 'Latex')
set(gca, 'FontSize', 15) 
axis tight
drawnow

%% PROOF 

if useintervals
    rho = intval(rho); % Careful, this is fine only because rho = 500 is an exact floating point number 
end

[rmin, rmax] = proof_SH(u, rho, nu, rstar); % This is where everything happens

if not(useintervals)
    fprintf("\nWARNING: This is not a fully rigorous proof, because you are not using interval arithmetic\n")
end



