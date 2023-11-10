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
rho=1025; % density of blood
rho=997; % density of blood

d=0.002; % pipe radius
U=0.5;

% % % Parameters dimensionless
Re=U*d*rho/mu;

% d=0.0015
% 
% u1=100*mu/(d*rho)
% u2=500*mu/(d*rho)

% 

% Re=400

U_scale= Re*mu/(d*rho);
gamma_scale= U_scale/d;

% 
% gamma_scale*3.7
% gamma_scale*1.3
% 
% gamma_scale*45

%

sr_vec=[logspace(0,5,2000)]';

x_current=load('data/fitted_parameters.txt');

% 
% %  Extract params from minima
p_cell=num2cell(x_current); % all of the free parameters
[alpha,beta,gamma_star,delta]=p_cell{:};
% parameters

xi=alpha*mu/(d^2*rho);
gamma_hat=gamma_star*d^2*rho/mu;
hat_beta=beta*mu/(d^2*rho);
LL=22.6;
delta;
% x_current=load('data/fit_params.txt')


% xi   = 0.0932
% kappa = 0.001 # diffusivity
% hatb =   1.2179e-04
% gamma_star= 2.3997e+04  #shear stress at which vwf unfolds
% delta=   7.8866e-06
% LL =       29.9996


x2=load('data/fitted_parameters.txt');


%  Extract params from minima
p_cell=num2cell(x2); % all of the free parameters
[alpha,beta,gamma_star,delta]=p_cell{:};


