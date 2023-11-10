%% This generates the all plots for main paper 
% This script generates all the plots for the main paper 


%---------------Assign plotting macros----------------------------
set(groot,'DefaultTextInterpreter','latex')
set(groot,'DefaultLegendInterpreter','latex')
set(groot,'DefaultAxesTickLabelInterpreter','latex')
set(groot,'DefaultTextFontname', 'CMU Serif')
set(groot,'DefaultAxesFontName', 'CMU Serif')
% set(groot,'DefaultAxesFontSize',11);
set(groot,'DefaultLineLineWidth',1)
set(0, 'DefaultFigureRenderer', 'painters');
set(groot,'DefaultAxesFontSize',12);
cmap_data = load('data/red_blue_cmap.mat');  %load in colourmap

parameters % import model parameters


%-------- Figure 4: VWF model in shear and elongational flow-------
% Note this figure is augmented in inkscape to add diagrams above.


if_mod=1;
fitted_params=load('data/fitted_parameters.txt'); % import fitted VWF parameters

sr_vec=logspace(0,5,1000); % array of shear rates to plot at 
params=[fitted_params,LL]; % add fixed VWF length  
[~,A]=vwf_extension_shear(params,sr_vec); %run VWF model shear flow
[~,A_e]=fene_extension_elong(params,if_mod,sr_vec); %run VWF model elongational flow

% Estimate 50% unfolding shear rate 
len=((A(:,1)+A(:,3))/2).^0.5; %calculate length 
half=max(((A(:,1)+A(:,3))/2).^0.5)/2; %calculate 50% of VWF's maximum length 

[C, ia, ~] = unique(len); % remove duplicate points where Len=1 or Len=LL
sr_half = interp1(C,sr_vec( ia),half); % find shear rate

% find values of extension for inset schematics
len_100 = interp1(sr_vec(ia),C,100);
len_2000 = interp1(sr_vec(ia),C,2e3);
len_5000 = interp1(sr_vec(ia),C,5e3);


len2=((A_e(:,1)+A_e(:,2))/2).^0.5; %calculate length 
[C, ia, ~] = unique(len2); 
sr_half_e = interp1(C,sr_vec(ia),half);

% find values of extension for inset schematics
len_100_e = interp1(sr_vec(ia),C,100);
len_2000_e = interp1(sr_vec(ia),C,2e3);
len_5000_e = interp1(sr_vec(ia),C,5e3);


% plot figure
figure1=figure('units','inch','position',[0,0,8.5,3.5]);
t=tiledlayout(1,2);
nexttile
plot(sr_vec,len,'r'); hold on
plot(sr_vec,(A(:,1)).^0.5,'k'); hold on
plot(sr_vec,(A(:,3)).^0.5,'color',[0.5,0.5,0.5])
set(gca,'XScale','log');
lk=xline(sr_half);
ylim([0,23])
scatter(100,len_100,'k','filled'); hold on
scatter(2e3,len_2000,'k','filled'); hold on
scatter(5e3,len_5000,'k','filled'); hold on


lk.LineStyle='-.';
legend({'$\mathcal{L}$','$\sqrt{A_{xx}}$','$\sqrt{A_{yy}}$','5,096s$^{-1}$'},'Location','northwest','FontSize',12);
ylabel('VWF Length' );
xlabel(' Shear rate $\dot{\gamma}$ (s$^{-1}$)');
xlim([100,1e5])
xticks([1,10^1,10^2,10^3,10^4,10^5])
if_mod=1;
ylim([0,23])

nexttile
plot(sr_vec,len2,'r'); hold on
plot(sr_vec,A_e(:,1).^0.5,'k','linewidth',1); hold on
plot(sr_vec,A_e(:,2).^0.5,'linewidth',1,'color',[0.5,0.5,0.5])
set(gca,'XScale','log');
xticks([1,10^1,10^2,10^3,10^4,10^5  ])
lk=xline(sr_half_e);
lk.LineStyle='-.';

scatter(100,len_100_e,'k','filled'); hold on
scatter(2e3,len_2000_e,'k','filled'); hold on
scatter(5e3,len_5000_e,'k','filled'); hold on

lk.Color='k';
legend({'$\mathcal{L}$','$\sqrt{A_{xx}}$',...
    '$\sqrt{A_{yy}}$','1,947s$^{-1}$'},'Location','northwest','FontSize',12);
xlabel(' Shear rate $\dot{\gamma}$ (s$^{-1}$)');
ylim([0,23])
xlim([100,1e5])

t.TileSpacing = 'compact';
t.Padding = 'compact';
box on
exportgraphics(figure1,'figs/fene.eps','ContentType','vector')

%%

%------Figure 6: heat plots of solution------------------------------------

% Import data for heat plots from selected simulation
Re=400;
h=0.5;
l1=1.5;
l2=2;
lc=0.025;

path='data/Re_'+string(Re)+'h-'+string(h)+'l1-'+string(l1)+'l2-'+string(l2)+'/';

 % Import all solutions as txt files
U=importdata(path+'usol.txt');
W=importdata(path+'wsol.txt');
E=importdata(path+'extsol.txt');
rr=importdata(path+'rrsol.txt');
zz=importdata(path+'zzsol.txt');
sr=importdata(path+'sr.txt');
psi=importdata(path+'psi.txt');
tot=importdata(path+'tot.txt');


% Triangulate data to enable plotting of meshed data 
bounda=0.5; %sensitivity for defining boundary from convex hull of mesh points. 

[~,~,~,~,zu,~]=extract(U,bounda);
[~,~,~,~,zw,~]=extract(W,bounda);
[~,~,~,~,ze,~]=extract(E,bounda);
[~,~,~,~,zrr,~]=extract(rr,bounda);
[DT,IO,x,y,zzz,k]=extract(zz,bounda);
[~,~,~,~,zsr,~]=extract(sr,bounda);
[~,~,~,~,zpsi,~]=extract(psi,bounda);
[~,~,~,~,ztot,~]=extract(tot,bounda);

%%Grids for getting streamlines
[xq,yq] = meshgrid(-10:0.01:30, 0:.001:1);
zpsi = griddata(x,y,zpsi,xq,yq);

% ze_cal=real(sqrt((zzz+zrr+1)/3)); 
% ze_cal=sqrt(((zzz+zrr)/2)+1)-1; % calculate the updated VWF extension definition
ze_cal=sqrt(ze+1)-1;% calculate the VWF extension 

painters_on=1; % ensures high qualitu figures, uncomment for speed.

heat_plot(x,y,(zu.^2+zw.^2).^1/2,xq,yq,zpsi,DT,IO,'(a) Fluid Velocity Magnitude','$|\bf{u}|$','U_plot',1,cmap_data,painters_on)
heat_plot(x,y,zsr,xq,yq,zpsi,DT,IO,'(b) Fluid Shear Rate','$\dot{\gamma}$','sr_plot',1,cmap_data,painters_on)
heat_plot(x,y,ze_cal,xq,yq,zpsi,DT,IO,'(c) VWF Extension','$\mathcal{E}$','E_plot',1,cmap_data,painters_on)


%% Extract data for flow types

Re=500;
h=0.5;
l1=1.5;
lc=0.025;

% case 1: l2=5 using strict boundary tolerance
l2=5;
path='data/Re_'+string(Re)+'h-'+string(h)+'l1-'+string(l1)+'l2-'+string(l2)+'/';
 tot=importdata(path+'tot.txt');
sr=importdata(path+'sr.txt');

rot=importdata(path+'rot.txt');
ex=importdata(path+'extsol.txt');

bounda=0.4;
[DT5,IO5,x5,y5,ztot5,k5]=extract(tot,bounda);
[~,~,~,~,sr5,~]=extract(sr,bounda);
[~,~,~,~,rot5,~]=extract(rot,bounda);

[~,~,~,~,ze5,~]=extract(ex,bounda);

% case 2: l2=2 using lenient boundary tolerance
 l2=2;
 path='data/Re_'+string(Re)+'h-'+string(h)+'l1-'+string(l1)+'l2-'+string(l2)+'/';
 tot2=importdata(path+'tot.txt');
 sr=importdata(path+'sr.txt');
 ex=importdata(path+'extsol.txt');
 rot=importdata(path+'rot.txt');

bounda=0.2;
[DT2,IO2,fitted_params,y2,ztot2,k2]=extract(tot2,bounda);
[~,~,~,~,sr2,~]=extract(sr,bounda);
[~,~,~,~,ze2,~]=extract(ex,bounda);
[~,~,~,~,rot2,~]=extract(rot,bounda);



%% Extract contour of 0.1 elongational flow

[xq,yq] = meshgrid(-10:0.1:10, -1:0.01:1);
vr = griddata([x5;x5],[y5;-y5],[ztot5;ztot5],xq,yq);
[c5]=contour(xq,yq,vr,-[0.2,2]); hold on

[xq,yq] = meshgrid(-10:0.1:10, -1:0.01:1);
vr = griddata([fitted_params;fitted_params],[y2;-y2],[ztot2;ztot2],xq,yq);
[c2]=contour(xq,yq,vr,-[0.2,-2]);

[start5,end5]=find(~c5(2,:));
 [start2,end2]=find(~c2);


%% Flow structure plots

cmap_data = load('red_blue_cmap.mat');  %load into structure

     set(0, 'DefaultFigureRenderer', 'opengl');
scale=1.3;
set(groot,'DefaultAxesFontSize',12*scale);

figure1=figure('units','inch','position',[0,0,8*scale,3*scale]);
t=tiledlayout(1,2);
nexttile % ax_vel = axes;

plot3(c5(1,end5(1):6:end5(2)),c5(2,end5(1):6:end5(2)),100*ones(size(c5(2,end5(1):6:end5(2)))),'k-.','LineWidth',0.2*scale); hold on

s=trisurf(DT5(IO5, :),x5,y5,(abs(sr5)-abs(rot5))./((abs(sr5)+abs(rot5))),'FaceColor','interp'); hold on
s.EdgeColor='none';
view(2)

% caxis([-1.2,1.2]);
 xlim([-10,10])
title('(a) Shallow Stenosis, $l_2=5$')
colormap(cmap_data.cmap)
plot3([-10,10,10,-10,-10],[0,0,1,1,0],[100,100,100,100,100],'k','LineWidth',0.2)
grid off
legend('Elongational flow: $(\dot{\gamma}-\dot{\omega})>0.1$')

 ylim([0,1])
    xlabel('$z$');
    ylabel('$r$','Interpreter','latex')

nexttile 
lam=(abs(sr2)-abs(rot2))./(abs(sr2)+abs(rot2));
plot3(c2(1,2:6:end),c2(2,2:6:end),100*ones(size(c2(1,2:6:end))),'k-.','LineWidth',0.2*scale); hold on
s=trisurf(DT2(IO2, :),fitted_params,y2,real(sqrt(lam)).*abs(sr2/2)*gamma_scale*1e-3,'FaceColor','interp'); hold on
s.EdgeColor='none';
view(2)
plot3([-10,10,10,-10,-10],[0,0,1,1,0],[100,100,100,100,100],'k','LineWidth',0.2)

%  caxis([-1.55,1.55]);
 ylim([0,1])
 xlim([-10,10])
 t.TileSpacing = 'compact';
t.Padding = 'compact';title('(b) Steep Stenosis, $l_2=2$')
a=colorbar('TickLabelInterpreter','latex');
a.Label.String = '$\dot{\gamma}-\dot{\omega}$';
a.Label.Interpreter = 'latex';
a.Label.FontSize = 1*scale;
box on
box on
    xlabel('$z$');
grid off
% exportgraphics(figure1,'data/class.eps','ContentType','vector')




cmap_data = load('red_blue_cmap.mat');  %load into structure

     set(0, 'DefaultFigureRenderer', 'opengl');
scale=1.3;
set(groot,'DefaultAxesFontSize',12*scale);

figure1=figure('units','inch','position',[0,0,8*scale,3*scale]);
t=tiledlayout(1,2);
nexttile % ax_vel = axes;

plot3(c5(1,end5(1):6:end5(2)),c5(2,end5(1):6:end5(2)),100*ones(size(c5(2,end5(1):6:end5(2)))),'k-.','LineWidth',0.2*scale); hold on

s=trisurf(DT2(IO2, :),fitted_params,y2,lam,'FaceColor','interp'); hold on
s.EdgeColor='none';
view(2)
colorbar
% caxis([-1.2,1.2]);
 xlim([-10,10])
title('(a) Shallow Stenosis, $l_2=5$')
colormap(cmap_data.cmap)
plot3([-10,10,10,-10,-10],[0,0,1,1,0],[100,100,100,100,100],'k','LineWidth',0.2)
grid off
legend('Elongational flow: $(\dot{\gamma}-\dot{\omega})>0.1$')

 ylim([0,1])
    xlabel('$z$');
    ylabel('$r$','Interpreter','latex')

nexttile 
lam=(abs(sr2)-abs(rot2))./(abs(sr2)+abs(rot2));
plot3(c2(1,2:6:end),c2(2,2:6:end),100*ones(size(c2(1,2:6:end))),'k-.','LineWidth',0.2*scale); hold on
s=trisurf(DT2(IO2, :),fitted_params,y2,real(sqrt(lam)).*abs(sr2/2)*gamma_scale*1e-3,'FaceColor','interp'); hold on
s.EdgeColor='none';
view(2)
plot3([-10,10,10,-10,-10],[0,0,1,1,0],[100,100,100,100,100],'k','LineWidth',0.2)

%  caxis([-1.55,1.55]);
 ylim([0,1])
 xlim([-10,10])
 t.TileSpacing = 'compact';
t.Padding = 'compact';title('(b) Steep Stenosis, $l_2=2$')
a=colorbar('TickLabelInterpreter','latex');
a.Label.String = '$\dot{\gamma}-\dot{\omega}$';
a.Label.Interpreter = 'latex';
a.Label.FontSize = 1*scale;
box on
box on
    xlabel('$z$');
grid off
% exportgraphics(figure1,'data/class.eps','ContentType','vector')




%% Find max SR and max E inside contour of elongation
 [in,~] = inpolygon(fitted_params,y2,c2(1,:),c2(2,:));

% max(sr2(in))
max(sr2)
% 
sqrt(max(ze2(in)))/21*100
sqrt(max(ze2))

[in,~] = inpolygon(x5,y5,c5(1,:),c5(2,:));

% max(sr5(in))
max(sr5)
sqrt(max(ze5(in)))/21*100
sqrt(max(ze5))/21*100

