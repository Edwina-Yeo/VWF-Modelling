close all

% minimise the distance of our model solution in shear flow to the Lippok
% et al data.
plot_display=0;

% [alpha,beta,gamma_star,delta,LL]
initial_guess=[];
initial_guess(1)=0.01;
initial_guess(2)=1e-3; % from Lippok et al.
initial_guess(3)=1e4; % from Lippok et al.
initial_guess(4)=1e-6; % 
LL=22.6;% from Schnieder et al. to obtain a maximum extension of 15, i.e. a max length of 16. 
% calculated: 16*sqrt(2)

lb=[1e-3,1e-8,4e3,1e-8]; % lower bound
ub=[0.1,1e-2,2e4,1e-3]; % upper bound
%
x2 = fmincon(@(params)vwf_fit_li(params,plot_display,LL),initial_guess,[],[],[],[],lb,ub);
%%
%  save('data/fitted_parameters.txt','x2','-ascii')

%

clear all
parameters
plot_display=1;

x_current=load('data/fitted_parameters.txt')


vwf_fit_li(x_current,plot_display,LL)



function error_sum=vwf_fit_li(params,plot_display,LL)


full_params=[params,LL];


sr_vec=[logspace(0,5,400)]';
lip=1./(1+exp(-(sr_vec-5522)/1271));

[length_v,~]=vwf_extension_shear(full_params,sr_vec);
rel_ext=(length_v-1)/max(length_v-1);%relative extension
rel_len=length_v/max(length_v);%relative extension


error=((rel_ext-lip));%./lip; % error at each point
error_sum=100*sum(abs(error),'all')./length(rel_ext) ;%percentage error
set(groot,'DefaultAxesFontSize',13);

if plot_display==1
figure1=figure('units','inch','position',[0,0,10,3.5]);
p_cell=num2cell(full_params); % all of the free parameters

[alpha,beta,gamma_star,delta,LL]=p_cell{:};

t=tiledlayout(1,3);
nexttile
 set(0, 'DefaultFigureRenderer', 'painters');
xlim([10,1e5])

 tau=alpha*((tanh(beta*(sr_vec-gamma_star))+1)/2+delta);

%  tau=alpha_art*((tanh(beta_art*(sr_vec-gamma_star_art))+1)/2+delta_art);
lk=xline(5096); hold on
lk.LineStyle='-.';
semilogx(sr_vec,tau,'k','linewidth',1); hold on
title('(a) Relaxation time $\tau$ (s)','Interpreter','latex')
set(gca,'XScale','log');
% xlim([100,5e4])
xticks([1,10,10^2,10^3,10^4,10^5])
ylim([0,0.08])

xlabel(' Shear rate $\dot{\gamma}$ (s$^{-1}$)','Interpreter','latex')
legend('$5,096$ s$^{-1}$','Location','northwest')
% exportgraphics(axes1,append('figs/relax.eps'),'Resolution',300) 

t.TileSpacing = 'compact';
t.Padding = 'compact';

box on
nexttile
        plot(sr_vec,rel_ext,'k'); hold on
plot(sr_vec,lip,'r');
   lk= xline(5096);
lk.LineStyle='-.';
xlim([10,1e5])

legend('VWF model','Lippok et al.','$5,096$ s$^{-1}$','Interpreter','latex','Location','northwest')
xlabel(' Shear rate $\dot{\gamma}$ (s$^{-1}$)');
ylim([0,1.1])
    set(gca,'XScale','log');
xticks([1,10,10^2,10^3,10^4,10^5])
title('(b) Relative VWF extension')


sch=[1, 1.0810810810810807;
10.00000000000001, 1.0810810810810807;
100, 1.1226611226611247;
398.1071705534969, 1.0810810810810807;
1751.7747448309788, 1.2058212058212057;
4863.585127384103, 12.515592515592514;
7864.156977348052, 15.01039501039501]; % extracted from Schnider extension data

nexttile
    plot(sr_vec,length_v,'k'); hold on
    plot(sch(:,1),sch(:,2),'r.','MarkerSize',10)

%              errorbar(exp_sr,(exp_ext),exp_vwf_err,'.','Color',  cmap(ceil(length(cmap)),:)); hold on
lk=xline(5096);
lk.LineStyle='-.';

xlim([10,1e5])
ylim([0,17])

legend('VWF model','Schneider et al.','$5,096$ s$^{-1}$','Interpreter','latex','Location','northwest')


% xlabel('Shear rate, ($s^{-1}$)')
% title('(b)')
title('(c) VWF length')

    set(gca,'XScale','log');
% exportgraphics(figure1 ,append('figs/model_vwf_fit.png'),'Resolution',300) 
xticks([1,10,10^2,10^3,10^4,10^5])
xlabel(' Shear rate $\dot{\gamma}$ ($s^{-1}$)');

exportgraphics(figure1,'figs/model_vwf_fit.eps','ContentType','vector')





figure
 set(0, 'DefaultFigureRenderer', 'painters');



sch=[1, 1.0810810810810807;
10.00000000000001, 1.0810810810810807;
100, 1.1226611226611247;
398.1071705534969, 1.0810810810810807;
1751.7747448309788, 1.2058212058212057;
4863.585127384103, 12.515592515592514;
7864.156977348052, 15.01039501039501]; % extracted from Schnider extension data

xlim([10,1e5])


        plot(sr_vec,rel_ext,'k'); hold on
plot(sr_vec,lip,'r');
   lk= xline(5096);
lk.LineStyle='-.';
xlim([10,1e5])

legend('VWF model','Lippok et. al','$5,096s^{-1}$','Interpreter','latex','Location','northwest')
xlabel(' Shear rate $\dot{\gamma}$ ($s^{-1}$)');
ylim([0,1.1])
    set(gca,'XScale','log');
xticks([1,10,10^2,10^3,10^4,10^5])
title('(b) Relative VWF extension')

nexttile
    plot(sr_vec,length_v,'k'); hold on
    plot(sch(:,1),sch(:,2),'r.','MarkerSize',10)

lk=xline(5096);
lk.LineStyle='-.';

xlim([10,1e5])
ylim([0,17])

legend('VWF model','Schneider et al.','$5,096s^{-1}$','Interpreter','latex','Location','northwest')

title('(c) VWF length')

    set(gca,'XScale','log');
xticks([1,10,10^2,10^3,10^4,10^5])
xlabel(' Shear rate $\dot{\gamma}$ ($s^{-1}$)');

exportgraphics(figure1,'figs/model_vwf_fit.png','ContentType','vector')
end

end

