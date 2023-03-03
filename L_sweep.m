%%
if_mod=1;

parameters

x2=load('data/fitted_parameters.txt')

set(groot,'DefaultLineLineWidth',1)
% vwf_params=[alpha,beta,gamma_star,delta,LL];
sr_vec=[logspace(0,5,2000)]';
if_mod=1

ll=LL./[7.5,5,2.5,1,1/2.5,1/5.5];
 figure1=figure('units','inch','position',[0,0,8.5,3.5]);

     set(0, 'DefaultFigureRenderer', 'opengl');

for i=1:length(ll)
params=[x2,ll(i)];
[len,A]=vwf_extension_shear(params,sr_vec); %run VWF model
[len2,A_e]=fene_extension_elong(params,if_mod,sr_vec); %run VWF model
half=max(((A(:,1)+A(:,3))/2).^0.5)/2; %21.0788/2
% Estimate middles
len=((A(:,1)+A(:,3))/2).^0.5;
[C, ia, ic] = unique(len);
sr_half = interp1(C,sr_vec(ia),half)
% lk=xline(sr_half);
half=max(((A_e(:,1)+A_e(:,2))/2).^0.5)/2; %21.0788/2

len2=((A_e(:,1)+A_e(:,2))/2).^0.5
[C, ia, ic] = unique(len2);
sr_half_e = interp1(C,sr_vec(ia),half)

% subplot(1,2,2)

%  Plot 2: Shear flow solutions in modified case-----------------------
set(groot,'DefaultAxesFontSize',12);
subplot(1,2,1)
plot(sr_vec,(((A(:,1)+A(:,3))/2).^0.5-1)./(max(((A(:,1)+A(:,3))/2).^0.5-1,[],'all'))); hold on

set(gca,'XScale','log');
legend({'$\mathcal{L}$','$\sqrt{A_{xx}}$','$\sqrt{A_{yy}}$','5,096s$^{-1}$'},'Location','northwest','FontSize',12);
ylabel('VWF Length' );
xlabel(' Shear rate $\dot{\gamma}$ ($s^{-1}$)');
xlim([100,1e5])
xticks([1,10^1,10^2,10^3,10^4,10^5 ...
    ])
if_mod=1;
set(groot,'DefaultLineLineWidth',1)

xlim([100,1e4])

% if i==1||i==length(ll)
% subplot(1,2,1)
% sr_half

    xline(sr_half)
%     i
%     20
% end

subplot(1,2,2)

plot(sr_vec,(((A_e(:,1)+A_e(:,2))/2).^0.5-1)./(max(((A_e(:,1)+A_e(:,2))/2).^0.5-1,[],'all'))); hold on
% plot(sr_vec_e,A_e(:,1).^0.5,'k','linewidth',1); hold on
% plot(sr_vec_e,zeros(size(A_e(:,1))),'k','linewidth',1); hold on

% plot(sr_vec_e,A_e(:,2).^0.5,'linewidth',1,'color',[0.5,0.5,0.5])

set(gca,'XScale','log');
% set(gca,'YScale','log');

% lk=xline(5e3);
% ylim([5e-4,1e4])
xticks([1,10^1,10^2,10^3,10^4,10^5 ...
    ])

% lk=xline(sr_half_e);
% lk.LineStyle='-.';
% lk.Color='b';
legend({'$\mathcal{L}$','$\sqrt{A_{xx}}$','$\sqrt{A_{yy}}$','1,947s$^{-1}$'},'Location','northwest','FontSize',12);

% legend({'$\mathcal{L}^2$','$A_{xx}$','$A_{yy}$','5,000s$^{-1}$'},'Location','northwest');
% title('Elongational flow' );
xlabel(' Shear rate $\dot{\gamma}$ ($s^{-1}$)');
% exportgraphics(axes1,append('figs/fene_elong.png'),'Resolution',300) 
% ylim([0,23])
xlim([100,1e4])

% if i==1||i==length(ll)
%     subplot(1,2,2)
% sr_half_e
xline(sr_half_e)
% i
% 10
% end
% t.TileSpacing = 'compact';
% t.Padding = 'compact';
box on
% exportgraphics(figure1,'figs/fene.eps','ContentType','vector')

end
