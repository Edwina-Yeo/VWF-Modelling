%% Parameters for Simulations 
% Author: Edwina Yeo 
% Contact: yeo@maths.ox.ac.uk 
% Updated 25/03/22


% Plotting parameters
set(groot,'DefaultTextInterpreter','latex')
set(groot,'DefaultLegendInterpreter','latex')
set(groot,'DefaultAxesTickLabelInterpreter','latex')
set(groot,'DefaultTextFontname', 'CMU Serif')
set(groot,'DefaultAxesFontName', 'CMU Serif')
set(groot,'DefaultAxesFontSize',11);
set(groot,'DefaultLineLineWidth',1)

% % Dimensional values
mu=    0.0025; %viscosity of blood
rho=997; % density of blood
d=0.002; % pipe radius
U=0.5;
% % Parameters dimensionless
Re=U*d*rho/mu
% 
xi=alpha*mu/(d^2*rho)
gamma_hat=gamma_star*d^2*rho/mu
hat_beta=beta*mu/(d^2*rho)
LL
 delta


