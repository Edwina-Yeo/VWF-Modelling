close all

% minimise the distance of our model solution in shear flow to the Lippok
% et al data.
plot_display=0;

% [alpha,beta,gamma_star,delta,LL]
initial_guess=[];
initial_guess(1)=1;
initial_guess(2)=2/1271; % from Lippok et al.
initial_guess(3)=5522; % from Lippok et al.
initial_guess(4)=1e-3; % 
initial_guess(5)=30; % from Schnieder et al.


lb=[0,0,1e3,0,10] % lower bound
ub=[10,1,1e6,0.1,100] % upper bound

% use multistart optimisation algorithm to find global minimas

parpool('local',2);
fun= @ (params)  vwf_fit_li(params,plot_display);
ms = MultiStart('UseParallel',true,'Display','final');
opts = optimoptions(@fmincon,'Algorithm','sqp');

problem = createOptimProblem('fmincon','objective',...
    fun,'x0',initial_guess,'lb',lb,'ub',ub,'options',opts); 

[xminm,fminm,flagm,outptm,manyminsm] = run(ms,problem,2000);  % run minimisation with this number of initial points

save(append('data/manyminsm.mat'),'manyminsm')  % Save the structure

delete(gcp('nocreate'))

%%
% 

parameters
load('data/manyminsm.mat')


fs=arrayfun(@(x)x.Fval,manyminsm);

%%

close all
plot(1:length(fs),fs)
fs(1)
plot_display=1

for j=20:30
    xmin=[];
for i=1:5
    xs=arrayfun(@(x)x.X(i),manyminsm);
    x_min(i)=xs(j);
end
x_min

vwf_fit_li(x_min,plot_display)
fs(j)

end


% x2 = fmincon(@(params)vwf_fit_li(params,plot_display),initial_guess,[],[],[],[],lb,ub)

%%

for i=1:5
    x2(i)
    initial_guess(i)
end

%%
x2=load('data/fit_params.txt')

%  Extract params from minima
p_cell=num2cell(x2); % all of the free parameters
[alpha,beta,gamma_star,delta,LL]=p_cell{:};
parameters

x2=load('data/fit_params.txt')


close all

plot_display=1;
vwf_fit_li(x2,plot_display)

% 
% 
% %%
% %  Extract params from minima
% p_cell=num2cell(x2); % all of the free parameters
% [alpha,beta,gamma_star,delta,LL]=p_cell{:};
% 
% % Plotting FENE P plots
% sr_vec=linspace(100,1e5,1e3); % shear rate vectors 
% set(groot,'DefaultAxesFontSize',10)
% set(groot,'DefaultLineLineWidth',1)
% 
% % VWF_params=[alpha,beta,gamma_star,delta,LL];
% figure1=figure('units','inch','position',[0,0,3,3]);
% axes1 = axes('Parent',figure1); 
% for i=1:1
% % Plot 1: Relaxation time dimensional ----------------------------------
%  tau=alpha*((tanh(beta*(sr_vec-gamma_star))+1)/2+delta);
% 
% %  tau=alpha_art*((tanh(beta_art*(sr_vec-gamma_star_art))+1)/2+delta_art);
% lk=xline(5e3); hold on
% lk.LineStyle='-.';
% semilogx(sr_vec,tau,'k','linewidth',1); hold on
% ylabel('Relaxation time $\tau$ (s)','Interpreter','latex')
% set(gca,'XScale','log');
% % xlim([100,5e4])
% 
% xlabel(' Shear rate $\dot{\gamma}$ ($s^{-1}$)','Interpreter','latex')
% legend('$5,000s^{-1}$','Location','northwest')
% exportgraphics(axes1,append('figs/relax.eps'),'Resolution',300) 
% end 
% 
% % 
% % 




function error_sum=vwf_fit_li(params,plot_display)

sr_vec=[logspace(0,3,50),logspace(3,5,200)]';
lip=1./(1+exp(-(sr_vec-5522)/1271));

[length_v,A]=vwf_extension_shear(params,sr_vec);
rel_ext=(length_v-1)/max(length_v-1);%relative extension
rel_len=length_v/max(length_v);%relative extension


length_v



error=((rel_ext-lip))./lip; % error at each point
error_sum=100*sum(abs(error),'all')./length(rel_ext) ;%percentage error
set(groot,'DefaultAxesFontSize',11);

if plot_display==1
figure1=figure('units','inch','position',[0,0,8,3]);
p_cell=num2cell(params); % all of the free parameters

[alpha,beta,gamma_star,delta,LL]=p_cell{:};

t=tiledlayout(1,3);
nexttile
 set(0, 'DefaultFigureRenderer', 'painters');
xlim([10,1e5])

 tau=alpha*((tanh(beta*(sr_vec-gamma_star))+1)/2+delta);

%  tau=alpha_art*((tanh(beta_art*(sr_vec-gamma_star_art))+1)/2+delta_art);
lk=xline(5148); hold on
lk.LineStyle='-.';
semilogx(sr_vec,tau,'k','linewidth',1); hold on
title('(a) Relaxation time $\tau$ (s)','Interpreter','latex')
set(gca,'XScale','log');
% xlim([100,5e4])
xticks([1,10,10^2,10^3,10^4,10^5])
ylim([0,0.17])

xlabel(' Shear rate $\dot{\gamma}$ ($s^{-1}$)','Interpreter','latex')
legend('$5,148s^{-1}$','Location','northwest')
% exportgraphics(axes1,append('figs/relax.eps'),'Resolution',300) 

t.TileSpacing = 'compact';
t.Padding = 'compact';

box on
nexttile
        plot(sr_vec,rel_ext,'k'); hold on
plot(sr_vec,lip,'r');
   lk= xline(5e3);
lk.LineStyle='-.';
xlim([10,1e5])

legend('VWF model','Lippok et. al','$5,148s^{-1}$','Interpreter','latex','Location','northwest')
xlabel(' Shear rate $\dot{\gamma}$ ($s^{-1}$)');
ylim([0,1.3])
    set(gca,'XScale','log');
xticks([1,10,10^2,10^3,10^4,10^5])
title('(b) Relative VWF extension')
nexttile
%   plot(sr_vec,lip); hold on
    plot(sr_vec,length_v,'k')
%              errorbar(exp_sr,(exp_ext),exp_vwf_err,'.','Color',  cmap(ceil(length(cmap)),:)); hold on
lk=xline(5e3);
lk.LineStyle='-.';

xlim([10,1e5])

legend('VWF model','$5,148s^{-1}$','Interpreter','latex','Location','northwest')

% xlabel('Shear rate, ($s^{-1}$)')
% title('(b)')
title('(c) VWF length')

    set(gca,'XScale','log');
% exportgraphics(figure1 ,append('figs/model_vwf_fit.png'),'Resolution',300) 
xticks([1,10,10^2,10^3,10^4,10^5])
xlabel(' Shear rate $\dot{\gamma}$ ($s^{-1}$)');

exportgraphics(figure1,'figs/model_vwf_fit.eps','ContentType','vector')
end

end

