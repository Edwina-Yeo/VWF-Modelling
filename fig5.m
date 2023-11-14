close all

% minimise the distance of our model solution in shear flow to the Lippok
% et al data.
plot_display=0;
LL=22.6;% from Schnieder et al. to obtain a maximum extension of 15, i.e. a max length of 16. 
% calculated: 16*sqrt(2)

lb=[1e-3,1e-8,4e3,1e-8]; % lower bound
ub=[0.1,1e-2,2e4,1e-3]; % upper bound
%

clc
plot_display=0

x_current=load('data/fitted_parameters.txt')
x2=x_current

ll=[5,10,15,50,75,100];
params_fitted_L=zeros(6,4)
for i=3:6
x2 = fmincon(@(params)vwf_fit_li(params,plot_display,ll(i)),x2,[],[],[],[],lb,ub)
params_fitted_L(i,:)=x2;
i
vwf_fit_li(x2,1,ll(i))

end

x2=x_current;
for i=3:-1:1
    x2 = fmincon(@(params)vwf_fit_li(params,plot_display,ll(i)),x2,[],[],[],[],lb,ub)
params_fitted_L(i,:)=x2;
i
vwf_fit_li(x2,1,ll(i))


end
% end
% 
% vwf_fit_li(x2,1,ll(i))
% 
% params=[x2,ll(i)]
% 
% subplot(1,3,2)
% sr_vec=[logspace(0,5,400)]';
% 
% [length_v,~]=vwf_extension_shear(params,sr_vec);
% rel_ext=(length_v-1)/max(length_v-1);%relative extension
% rel_len=length_v/max(length_v);%relative extension
% 
%         plot(sr_vec,rel_ext,'b'); hold on



%%
close all
ll=[5,10,15,50,75,100];
parameters

figure1=figure('units','inch','position',[0,0,6,3]);
% save("params_fitted_L")

params_fitted=load('params_fitted_L' )
params_fitted_L=params_fitted.params_fitted_L
% if_mod=1;
% for i=1:6
%     params=[params_fitted_L(i,:),ll(i)]
%     vwf_fit_li(params,1,ll(i))
% sr_vec=[logspace(0,5,2000)]';
% [len2,A_e]=fene_extension_elong(params,if_mod,sr_vec); %run VWF model
% rel_ext_e=(len2-1)/max(len2-1);%relative extension
% [length_v,~]=vwf_extension_shear(params,sr_vec);
% rel_ext=(length_v-1)/max(length_v-1);%relative extension
% rel_len=length_v/max(length_v);%relative extension
% subplot(1,2,1)
% lip=1./(1+exp(-(sr_vec-5522)/1271));
% plot(sr_vec,lip,'r')
%         C=plot(sr_vec,rel_ext); hold on
% %         clabel(C,'l=')
%                 xlim([100,1e5])
% 
% legend('VWF model','Lippok et. al','5,096 s$^{-1}$','Interpreter','latex','Location','northwest')
% xlabel(' Shear rate $\dot{\gamma}$ ($s^{-1}$)');
% ylim([0,1.1])
%     set(gca,'XScale','log');
% xticks([1,10,10^2,10^3,10^4,10^5])
% subplot(1,2,2)
% plot(sr_vec,rel_ext_e); hold on
%         xlim([100,1e5])
% 
% legend('VWF model','Lippok et. al','$5,096s^{-1}$','Interpreter','latex','Location','northwest')
% xlabel(' Shear rate $\dot{\gamma}$ ($s^{-1}$)');
% ylim([0,1.1])
%     set(gca,'XScale','log');
% xticks([1,10,10^2,10^3,10^4,10^5])
% end


%
close all

i=1

params=[params_fitted_L(i,:),ll(i)]
sr_vec=[logspace(0,5,2000)]';
[len2,A_e]=fene_extension_elong(params,if_mod,sr_vec); %run VWF model
rel_ext_e1=(len2-1)/max(len2-1);%relative extension
[length_v,~]=vwf_extension_shear(params,sr_vec);
rel_ext1=(length_v-1)/max(length_v-1);%r
% 
% 
% 

i=6
params=[params_fitted_L(i,:),ll(i)]
sr_vec=[logspace(0,5,2000)]';
[len2,A_e]=fene_extension_elong(params,if_mod,sr_vec); %run VWF model
rel_ext_e6=(len2-1)/max(len2-1);%relative extension
[length_v,~]=vwf_extension_shear(params,sr_vec);
rel_ext6=(length_v-1)/max(length_v-1);%relative extension



%

half=0.5; %21.0788/2


[C, ia, ic] = unique(rel_ext_e6);
sr_half6 = interp1(C,sr_vec(ia),half);
lk6=xline(sr_half6);


[C, ia, ic] = unique(rel_ext_e1);
sr_half1 = interp1(C,sr_vec(ia),half);
lk1=xline(sr_half1);
%
figure
close all
clc


params=[x_current,LL];
[len2,A_e]=fene_extension_elong(params,if_mod,sr_vec); %run VWF model
rel_ext_e=(len2-1)/max(len2-1);%relative extension
[length_v,~]=vwf_extension_shear(params,sr_vec);
rel_ext=(length_v-1)/max(length_v-1);%r
%
 set(0, 'DefaultFigureRenderer', 'painters');


close all
figure1=figure('units','inch','position',[0,0,8.25,4 ...
    ]);
t=tiledlayout(1,2)
x2 = [sr_vec; flipud(sr_vec)];
inBetweenshear = [rel_ext1; flipud(rel_ext6)];
inBetween_elong = [rel_ext_e1;flip(rel_ext_e6)];

nexttile

lip=1./(1+exp(-(sr_vec-5522)/1271));
l=patch(x2, inBetweenshear, [0.8,0.8,0.8]); hold on
l.EdgeColor='none';
plot(sr_vec,lip,'r'); hold on
        C=plot(sr_vec,rel_ext,'k'); hold on
%         clabel(C,'l=')
                xlim([200,5e4])
                lk=xline(5096);
lk.LineStyle='-.';
box on 
title('(a) Fitted pure shear flow behaviour')

legend('Range of fit for L = 5 - 100','Lippok et al.','Model fit for L $=22.6$','5,096s$^{-1}$','Interpreter','latex','Location','northwest')
xlabel('Shear rate $\dot{\gamma}$ (s$^{-1}$)');
ylabel('VWF extension, $\mathcal{E}$');

ylim([0,1.1])
    set(gca,'XScale','log');
xticks([1,10,10^2,10^3,10^4,10^5])
nexttile
l=patch(x2, inBetween_elong, [0.8,0.8,0.8]); hold on
l.EdgeColor='none';
plot(sr_vec,rel_ext_e,'k'); hold on
                xlim([200,5e4])

lk1=xline(sr_half6);
lk1.LineStyle='-.';
lk1=xline(sr_half1);
lk1.LineStyle='--';
lk=legend('Prediction range','Prediction for L $=22.6$','Min $=635$ s$^{-1}$ (L $=100$)','Max $=3,280$ s$^{-1}$ (L $=5$)','Interpreter','latex','Location','northwest')
xlabel(' Shear rate $\dot{\gamma}$ (s$^{-1}$)');
% ylabel('VWF extension, $\mathcal{E}$)');

lk.Location='southeast';
title('(b) Predicted pure elongational flow behaviour')
ylim([0,1.1])
    set(gca,'XScale','log');
xticks([1,10,10^2,10^3,10^4,10^5])


box on

t.TileSpacing='compact';
t.Padding='tight';
exportgraphics(figure1,'figs/model_params.eps','ContentType','vector')


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
% figure1=figure('units','inch','position',[0,0,10,3.5]);
p_cell=num2cell(full_params); % all of the free parameters

[alpha,beta,gamma_star,delta,LL]=p_cell{:};

% t=tiledlayout(1,3)
subplot(1,3,1)
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

xlabel(' Shear rate $\dot{\gamma}$ ($s^{-1}$)','Interpreter','latex')
legend('$5,096s^{-1}$','Location','northwest')
% exportgraphics(axes1,append('figs/relax.eps'),'Resolution',300) 

subplot(1,3,2)

        plot(sr_vec,rel_ext,'k'); hold on
plot(sr_vec,lip,'r');
   lk= xline(5096);
lk.LineStyle='-.';
xlim([10,1e5])

legend('VWF model','Lippok et al.','$5,096s^{-1}$','Interpreter','latex','Location','northwest')
xlabel(' Shear rate $\dot{\gamma}$ ($s^{-1}$)');
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

subplot(1,3,3)

    plot(sr_vec,length_v,'k'); hold on
    plot(sch(:,1),sch(:,2),'r.','MarkerSize',10)

%              errorbar(exp_sr,(exp_ext),exp_vwf_err,'.','Color',  cmap(ceil(length(cmap)),:)); hold on
lk=xline(5096);
lk.LineStyle='-.';

xlim([10,1e5])
ylim([0,17])

legend('VWF model','Schneider et al.','$5,096s^{-1}$','Interpreter','latex','Location','northwest')


% xlabel('Shear rate, ($s^{-1}$)')
% title('(b)')
title('(c) VWF length')

    set(gca,'XScale','log');
% exportgraphics(figure1 ,append('figs/model_vwf_fit.png'),'Resolution',300) 
xticks([1,10,10^2,10^3,10^4,10^5])
xlabel(' Shear rate $\dot{\gamma}$ ($s^{-1}$)');


end

end

