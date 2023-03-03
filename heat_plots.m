% 


set(groot,'DefaultTextInterpreter','latex')
set(groot,'DefaultLegendInterpreter','latex')
set(groot,'DefaultAxesTickLabelInterpreter','latex')
set(groot,'DefaultTextFontname', 'CMU Serif')
set(groot,'DefaultAxesFontName', 'CMU Serif')
set(groot,'DefaultAxesFontSize',11);
set(groot,'DefaultLineLineWidth',1)


%% 
% run with re=500, h=0.5, l2=2, l1=2, l=30+l1+l2

maxe=load('data/Convergence/max_e.txt')
mine=load('data/Convergence/min_e.txt')
maxw=load('data/Convergence/max_w.txt')
verts=load('data/Convergence/meshverts.txt')

maxe2=load('data/Convergence/max_e2.txt')
mine2=load('data/Convergence/min_e2.txt')
maxw2=load('data/Convergence/max_w2.txt')
verts2=load('data/Convergence/meshverts2.txt')
maxe3=load('data/Convergence/max_e3.txt')
mine3=load('data/Convergence/min_e3.txt')
maxw3=load('data/Convergence/max_w3.txt')
verts3=load('data/Convergence/meshverts3.txt')

maxe4=load('data/Convergence/max_e4.txt')
mine4=load('data/Convergence/min_e4.txt')
maxw4=load('data/Convergence/max_w4.txt')
verts4=load('data/Convergence/meshverts4.txt')

% maxe=maxe(3:end)
% mine=mine(2:end)
% maxw=maxw(2:end)
mine=[mine(2:3)',mine4,mine2,mine3,mine(4:end)']'
maxe=[maxe(2:3)',maxe4,maxe2,maxe3,maxe(4:end)']'
maxw=[maxw(2:3)',maxw4,maxw2,maxw3,maxw(4:end)']'
verts=[verts(2:3)',verts4,verts2,verts3,verts(4:end)']'


%%


steps=[100,250,750,1000,1500]

min_step=[]
max_step=[]
for i=1:length(steps)
    i
min_step(i)=load(append('data/Convergence/Re_500h-0.5l1-1.5l2-1.5steps',string(steps(i)),'min_e_step.txt'))
max_step(i)=load(append('data/Convergence/Re_500h-0.5l1-1.5l2-1.5steps',string(steps(i)),'max_e_step.txt'))
end


%%
% figure
% loglog(verts(3:end),abs((maxe(3:end)-maxe(2)))/(maxe(2)),'.-'); hold on
% loglog(verts(3:end),abs((maxw(3:end)-maxw(2)))/(maxw(2)),'.-')
% % loglog(verts(3:end),abs((mine(3:end))),'.-')

set(groot,'DefaultAxesFontSize',11);

lcs=[0.4,0.3,0.2,0.1,0.07,0.05,0.02,0.01]


close all
loc=4

Nplot=logspace(1,5,10);

figure1=figure('units','inch','position',[0,0,6,3]);
t=tiledlayout(1,2)
nexttile
% loglog(times,abs(maxs-max_pb),'k.-','MarkerSize',10) ;hold on
 loglog(verts(2:end),abs((maxe(2:end)-maxe(1)))/(maxe(1)),'k.-','MarkerSize',10) ;hold on
loglog(verts(2:end),abs((maxw(2:end)-maxw(1)))/(maxw(1)),'.-','MarkerSize',10,'Color',[0.5,0.5,0.5]) ;hold on
 loglog(verts(loc),abs((maxe(loc)-maxe(1)))/(maxe(1)),'r.-','MarkerSize',10) ;hold on
loglog(verts(loc),abs((maxw(loc)-maxw(1)))/(maxw(1)),'r.-','MarkerSize',10) ;hold on
Nplot=logspace(3.8,5,10);
%  loglog(Nplot,(20*Nplot.^(-1)),'k-.','MarkerSize',10); hold on

Nplot=logspace(3.8,5,10);
%  loglog(Nplot,(5000*Nplot.^(-2)),'k-.','MarkerSize',10,'Color',[0.5,0.5,0.5]) ;hold on
%  loglog(Nplot,(Nplot.^(-1)),'b-','MarkerSize',10); hold on

% loglog(times(4),abs(maxs(4)-max_pb),'r.','MarkerSize',10) ;hold on

xlabel('Mesh vertices')
ylabel('Relative error in maximum')

title('(a) Maximum convergence')
xticks([1e2,1e3,1e4,1e5])
yticks([1e-5,1e-4,1e-3,1e-2,1e-1])

legend('VWF Eq.','Navier Stokes Eq.','Selected mesh','Location','southwest')
% exportgraphics(figure1,'data/validation/time_validation.png','Resolution',300) 

% ylim([10^(-3),20])
% xlim([10^(-3.8),0.3])
% exportgraphics(figure1,'data/mesh_validation.eps','ContentType','vector')



% %%
% close all

nexttile
% Nplot=logspace(-2.8,-1.8,10);

% figure1=figure('units','inch','position',[0,0,3,3]);
% loglog(times,abs(maxs-max_pb),'k.-','MarkerSize',10) ;hold on
%  loglog(Nplot,(2*Nplot),'b-','MarkerSize',10)
 semilogx(verts(2:end),mine(2:end)/maxe(1),'k.-','MarkerSize',10) ;hold on
 semilogx(verts(loc),mine(loc)/maxe(1),'r.-','MarkerSize',10) ;hold on

%  loglog(verts(end-2),abs((maxe(end-2)-maxe(2)))/(maxe(2)),'r.-','MarkerSize',10) ;hold on

% % loglog(times(4),abs(maxs(4)-max_pb),'r.','MarkerSize',10) ;hold on
% yline(1)
xlabel('Mesh vertices')
ylabel('Minimum extension (normalised)')
xticks([1e2,1e3,1e4,1e5])
% yticks([1e-5,1e-4,1e-3,1e-2,1e-1])
title('(b) Minimum extension convergence')

legend('VWF Eq.','Selected mesh','Location','southeast')
% exportgraphics(figure1,'data/validation/time_validation.png','Resolution',300) 
% yticks([1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2])

% ylim([10^(-3),20])
% xlim([10^(-3.8),0.3])
% exportgraphics(figure1,'data/min_fene.png','Resolution',300) 

t.TileSpacing = 'compact';
t.Padding = 'compact';
exportgraphics(figure1,'data/convergence.eps','ContentType','vector')

%%
close all


% Nplot=logspace(-2.8,-1.8,10);

figure1=figure('units','inch','position',[0,0,3,3]);
% loglog(times,abs(maxs-max_pb),'k.-','MarkerSize',10) ;hold on
%  loglog(Nplot,(2*Nplot),'b-','MarkerSize',10)
 loglog(verts(2:end),100*abs(mine(2:end))/maxe(2),'k.-','MarkerSize',10) ;hold on
%   loglog(verts(loc),abs(mine(loc))/maxe(2),'r.-','MarkerSize',10) ;hold on

%  loglog(verts(end-2),abs((maxe(end-2)-maxe(2)))/(maxe(2)),'r.-','MarkerSize',10) ;hold on

% loglog(times(4),abs(maxs(4)-max_pb),'r.','MarkerSize',10) ;hold on
yline(1)
title('(b) Minimum convergence')

xlabel('Mesh vertices')
ylabel('Minimum extension (normalised)')
xticks([1e2,1e3,1e4,1e5])
% yticks([1e-5,1e-4,1e-3,1e-2,1e-1])

legend('VWF Eq.','Selected mesh','Location','southwest')
% exportgraphics(figure1,'data/validation/time_validation.png','Resolution',300) 
yticks([1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2])

% ylim([10^(-3),20])


t.TileSpacing = 'compact';
t.Padding = 'compact';
% xlim([10^(-3.8),0.3])
exportgraphics(figure1,'data/min_fene.eps','ContentType','vector')

% exportgraphics(figure1,'data/min_fene.png','Resolution',300) 
%%


parameters

x2=load('data/fitted_parameters.txt')

%  Extract params from minima
p_cell=num2cell(x2); % all of the free parameters
[alpha,beta,gamma_star,delta]=p_cell{:};


%%
if_mod=1;

parameters
x2=load('data/fitted_parameters.txt')

set(groot,'DefaultLineLineWidth',1)
% vwf_params=[alpha,beta,gamma_star,delta,LL];
sr_vec=[logspace(0,3,50),logspace(3,5,200)]';
if_mod=1

params=[x2,LL]
[len,A]=vwf_extension_shear(params,sr_vec) %run VWF model


[len,A_e,sr_vec_e]=fene_extension_elong(params,if_mod); %run VWF model


%% Estimate middles
len=((A(:,1)+A(:,3))/2).^0.5;
[C, ia, ic] = unique(len);
sr_half = interp1(C,sr_vec(ia),half)


len2=((A_e(:,1)+A_e(:,2))/2).^0.5;
[C, ia, ic] = unique(len2);
sr_half_e = interp1(C,sr_vec_e(ia),half)

%%



%  Plot 2: Shear flow solutions in modified case-----------------------
close all
half=max(((A(:,1)+A(:,3))/2).^0.5)/2 %21.0788/2
 set(0, 'DefaultFigureRenderer', 'painters');
set(groot,'DefaultAxesFontSize',12);
 figure1=figure('units','inch','position',[0,0,8.5,3.5]);
t=tiledlayout(1,2);
nexttile
plot(sr_vec,((A(:,1)+A(:,3))/2).^0.5,'r'); hold on
plot(sr_vec,(A(:,1)).^0.5,'k'); hold on
% plot(sr_vec,A(:,2),'k-.')
plot(sr_vec,(A(:,3)).^0.5,'color',[0.5,0.5,0.5])
set(gca,'XScale','log');
% set(gca,'YScale','log');
lk=xline(sr_half);
% ylim([5e-4,1e4])
 ylim([0,23])

% xlim([100,5e4])
% xticks([10^2,10^3,10^4])

lk.LineStyle='-.';
legend({'$\mathcal{L}$','$\sqrt{A_{xx}}$','$\sqrt{A_{yy}}$','5,096s$^{-1}$'},'Location','northwest','FontSize',12);
ylabel('VWF Length' );
xlabel(' Shear rate $\dot{\gamma}$ ($s^{-1}$)');
% title('Shear flow' );
xlim([100,1e5])
% ylim([5e-4,1e4])
xticks([1,10^1,10^2,10^3,10^4,10^5 ...
    ])
if_mod=1;
set(groot,'DefaultLineLineWidth',1)
ylim([0,23])
% xlim([100,1e4])

nexttile
xlim([100,1e4])

plot(sr_vec_e,((A_e(:,1)+A_e(:,2))/2).^0.5,'r'); hold on
plot(sr_vec_e,A_e(:,1).^0.5,'k','linewidth',1); hold on
% plot(sr_vec_e,zeros(size(A_e(:,1))),'k','linewidth',1); hold on

plot(sr_vec_e,A_e(:,2).^0.5,'linewidth',1,'color',[0.5,0.5,0.5])

set(gca,'XScale','log');
% set(gca,'YScale','log');

% lk=xline(5e3);
% ylim([5e-4,1e4])
xticks([1,10^1,10^2,10^3,10^4,10^5 ...
    ])

lk=xline(sr_half_e);
lk.LineStyle='-.';
lk.Color='b';
legend({'$\mathcal{L}$','$\sqrt{A_{xx}}$','$\sqrt{A_{yy}}$','1,947s$^{-1}$'},'Location','northwest','FontSize',12);

% legend({'$\mathcal{L}^2$','$A_{xx}$','$A_{yy}$','5,000s$^{-1}$'},'Location','northwest');
% title('Elongational flow' );
xlabel(' Shear rate $\dot{\gamma}$ ($s^{-1}$)');
% exportgraphics(axes1,append('figs/fene_elong.png'),'Resolution',300) 
ylim([0,23])
xlim([100,1e5])

t.TileSpacing = 'compact';
t.Padding = 'compact';
box on
exportgraphics(figure1,'figs/fene.eps','ContentType','vector')



%% Import data for heat plots
Re=500
h=0.5
l1=1.5
l2=2
lc=0.025

tdata = load('red_blue_cmap.mat');  %load into structure
path='data/Re_'+string(Re)+'h-'+string(h)+'l1-'+string(l1)+'l2-'+string(l2)+'/'
% path='data/Re_'+string(Re)+'h-'+string(h)+'l1-'+string(l1)+'l2-'+string(l2)+'lc'+string(lc)+'/'

 % Import all solutions as txt files
  U=importdata(path+'usol.txt');
 W=importdata(path+'wsol.txt');
 E=importdata(path+'extsol.txt');
 rr=importdata(path+'rrsol.txt');
 zz=importdata(path+'zzsol.txt');

 rot=importdata(path+'rot.txt');
 sr=importdata(path+'sr.txt');
 psi=importdata(path+'psi.txt');
 tot=importdata(path+'tot.txt');
%% Import data from FENICS solution and triangulate it
 bounda=0.5

[DT,IO,x,y,zu,k]=extract(U,bounda);
[~,~,~,~,zw,~]=extract(W,bounda);
[DT,IO,x,y,ze,k]=extract(E,bounda);

[DT,IO,x,y,zrr,k]=extract(rr,bounda);
[DT,IO,x,y,zzz,k]=extract(zz,bounda);
%%
[~,~,~,~,zsr,~]=extract(sr,bounda);
[~,~,~,~,zpsi,~]=extract(psi,bounda);
[~,~,~,~,ztot,~]=extract(tot,bounda);

%%Grids for getting streamlines
[xq,yq] = meshgrid(-10:0.01:30, 0:.001:1);
zpsi = griddata(x,y,zpsi,xq,yq);

%%

ze_cal=real(sqrt((zzz+zrr+1)/3))
painters_on=0;
heat_plot(x,y,(zu.^2+zw.^2).^1/2,xq,yq,zpsi,DT,IO,'(a) Fluid Velocity Magnitude','$|\bf{u}|$','U_plot',1,tdata,painters_on)
heat_plot(x,y,ze_cal,xq,yq,zpsi,DT,IO,'(c) VWF Extension','$\mathcal{E}$','E_plot',1,tdata,painters_on)
heat_plot(x,y,zsr,xq,yq,zpsi,DT,IO,'(b) Fluid Shear Rate','$\dot{\gamma}$','sr_plot',1,tdata,painters_on)

% %% Grids for getting orientation
% [xq,yq] = meshgrid(-10:0.5:30, 0.05:0.1:0.95);
% vr = griddata(x,y,zr,xq,yq);
% vz = griddata(x,y,zz,xq,yq);
% polyin = polyshape(x(k),y(k));
% 
% xflat=reshape(xq,1,[]);
% yflat=reshape(yq,1,[]);
% [TFin,TFon] = isinterior(polyin,xflat,yflat);
% indIN=reshape(TFin,size(xq));
% indon=reshape(TFon,size(xq));
% indon=not(indon);
% 
% AR=vr.*indIN.*indon-1;
% AR(AR==-1)=0;
% AR(AR==nan)=0;
% AZ=vz.*indIN.*indon-1;
% AZ(AZ==-1)=0;
% AZ(AZ==nan)=0;
% 
% set(0, 'DefaultFigureRenderer', 'painters');
% 
% 
% figure1=figure('units','inch','position',[0,0,7,2]);
% ax_vel = axes;
% s=trisurf(DT(IO, :),x,y,-ztot,'FaceColor','interp'); hold on
% s.EdgeColor='none';
% view(2)
% % % % % xlim([-10,30]);
% ylabel(ax_vel,'$r$','Interpreter','latex')
% a=colorbar(ax_vel,'TickLabelInterpreter','latex');
% a.Label.String = ' $\dot{\gamma}-\dot{\omega}$';
% a.Label.Interpreter = 'latex';
% % a.Limits=[-1.92,1.92];
% caxis([-1.92,1.92]);
% 
% colormap(tdata.cmap)
% plot3([-10,30,30,-10],[0,0,1,1],[100,100,100,100],'k','LineWidth',0.2)
% 
% min(tot)
% max(tot)
% title(ax_vel,'Flow classification','Interpreter','latex')
% % qH2=quiver3(xq,yq,50*ones(size(xq)),sqrt(abs(AZ)),-sqrt(abs(AR)),2*ones(size(xq)),5); hold on
% %  qH2.ShowArrowHead = 'off';
% %  qH2.Color = 'k';
%  ylim([0,1])
%  xlim([-10,30])
% 
%  xlabel('$z$','Interpreter','latex')
% ylabel('$r$','Interpreter','latex') 
% grid off
% 
% % colormap(tdata.cmap) 
% % exportgraphics(figure1,'data/orient4.png','Resolution',300) 
% plot3([-10,30,30,-10,-10],[0,0,1,1,0],[100,100,100,100,100],'k','LineWidth',0.2)
% box on
% % exportgraphics(figure1,'data/orient2.png','Resolution',300) 


%% Extract data for flow types

Re=500
h=0.5
l1=1.5
lc=0.025

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
[DT2,IO2,x2,y2,ztot2,k2]=extract(tot2,bounda);
[~,~,~,~,sr2,~]=extract(sr,bounda);
[~,~,~,~,ze2,~]=extract(ex,bounda);
[~,~,~,~,rot2,~]=extract(rot,bounda);



%% Extract contour of 0.1 elongational flow

[xq,yq] = meshgrid(-10:0.1:10, -1:0.01:1);
vr = griddata([x5;x5],[y5;-y5],[ztot5;ztot5],xq,yq);
[c5]=contour(xq,yq,vr,-[0.2,2]); hold on

[xq,yq] = meshgrid(-10:0.1:10, -1:0.01:1);
vr = griddata([x2;x2],[y2;-y2],[ztot2;ztot2],xq,yq);
[c2]=contour(xq,yq,vr,-[0.2,-2]);

[start5,end5]=find(~c5(2,:));
 [start2,end2]=find(~c2);


%% Flow structure plots

tdata = load('red_blue_cmap.mat');  %load into structure

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
colormap(tdata.cmap)
plot3([-10,10,10,-10,-10],[0,0,1,1,0],[100,100,100,100,100],'k','LineWidth',0.2)
grid off
legend('Elongational flow: $(\dot{\gamma}-\dot{\omega})>0.1$')

 ylim([0,1])
    xlabel('$z$');
    ylabel('$r$','Interpreter','latex')

nexttile 
lam=(abs(sr2)-abs(rot2))./(abs(sr2)+abs(rot2));
plot3(c2(1,2:6:end),c2(2,2:6:end),100*ones(size(c2(1,2:6:end))),'k-.','LineWidth',0.2*scale); hold on
s=trisurf(DT2(IO2, :),x2,y2,real(sqrt(lam)).*abs(sr2/2)*gamma_scale*1e-3,'FaceColor','interp'); hold on
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




tdata = load('red_blue_cmap.mat');  %load into structure

     set(0, 'DefaultFigureRenderer', 'opengl');
scale=1.3;
set(groot,'DefaultAxesFontSize',12*scale);

figure1=figure('units','inch','position',[0,0,8*scale,3*scale]);
t=tiledlayout(1,2);
nexttile % ax_vel = axes;

plot3(c5(1,end5(1):6:end5(2)),c5(2,end5(1):6:end5(2)),100*ones(size(c5(2,end5(1):6:end5(2)))),'k-.','LineWidth',0.2*scale); hold on

s=trisurf(DT2(IO2, :),x2,y2,lam,'FaceColor','interp'); hold on
s.EdgeColor='none';
view(2)
colorbar
% caxis([-1.2,1.2]);
 xlim([-10,10])
title('(a) Shallow Stenosis, $l_2=5$')
colormap(tdata.cmap)
plot3([-10,10,10,-10,-10],[0,0,1,1,0],[100,100,100,100,100],'k','LineWidth',0.2)
grid off
legend('Elongational flow: $(\dot{\gamma}-\dot{\omega})>0.1$')

 ylim([0,1])
    xlabel('$z$');
    ylabel('$r$','Interpreter','latex')

nexttile 
lam=(abs(sr2)-abs(rot2))./(abs(sr2)+abs(rot2));
plot3(c2(1,2:6:end),c2(2,2:6:end),100*ones(size(c2(1,2:6:end))),'k-.','LineWidth',0.2*scale); hold on
s=trisurf(DT2(IO2, :),x2,y2,real(sqrt(lam)).*abs(sr2/2)*gamma_scale*1e-3,'FaceColor','interp'); hold on
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
 [in,~] = inpolygon(x2,y2,c2(1,:),c2(2,:));

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
function heat_plot(x,y,zheat,xq,yq,zpsi,DT,IO,title_label,col_label,save_name,z_label_on,tdata,painters_on)
if painters_on
 set(0, 'DefaultFigureRenderer', 'painters');
else
    set(0, 'DefaultFigureRenderer', 'opengl');
end
set(groot,'DefaultAxesFontSize',12);

figure1=figure('units','inch','position',[0,0,8,2]);
ax_vel = axes;
ax_vel.Box='on';
ax_vel.YGrid='off';
ax_vel.XGrid='off';
ax_vel.ZGrid='off';
title(ax_vel,title_label,'Interpreter','latex')

s=trisurf(DT(IO, :),x,y,zheat,'FaceColor','interp'); hold on
s.EdgeColor='none';
view(2)
xlim([-10,30]);
if z_label_on
    xlabel('$z$');
end
plot3([-10,30,30,-10,-10],[0,0,1,1,0],[100,100,100,100,100],'k','LineWidth',0.2)
ylabel(ax_vel,'$r$','Interpreter','latex')
a=colorbar(ax_vel,'TickLabelInterpreter','latex');
a.Label.String = col_label;
a.Label.Interpreter = 'latex';
colormap(tdata.cmap)
title(ax_vel,title_label,'Interpreter','latex')
ax_vel2 = axes('Position',[ax_vel.Position],'Xlim',[ax_vel.XLim],'Ylim',[ax_vel.YLim],'PlotBoxAspectRatio',[ax_vel.PlotBoxAspectRatio]);
contour(ax_vel2,xq,yq,zpsi+10,10+[-0.1:0.01:-0.08,-0.08:0.02:0.25],'k'); hold on

linkaxes([ax_vel,ax_vel2]);
ax_vel2.Visible = 'off';
drawnow
grid off
ax_vel2.Box='on';

ax_vel.YGrid='off';
ax_vel.XGrid='off';
ax_vel.ZGrid='off';
ax_vel.TickDir='out';

title(ax_vel2,title_label,'Interpreter','latex')
exportgraphics(figure1,append('figs/',save_name,'.eps'),'ContentType','vector')


end

 %Extraction and plotting function for heatmaps of concentration:
%triangulation is required to avoid plotting over the aggregate
function [DT,IO,x,y,z,k]=extract(c,bounda)
%import data both c surface and quiver plots
%extract variables
x=c(:,1);
y=c(:,2);
z=c(:,3);
%extract aggregate boundary;
 k=boundary(x,y,bounda);
 C=[k(1:end-1),k(2:end)]; %form constraint for plotting
 DT=delaunayTriangulation([x,y],C); %create triangulation for plotting
 IO = isInterior(DT); %only plot inside the boundary.
end
