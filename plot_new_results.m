

set(groot,'DefaultAxesFontSize',12);

close all
figure1=figure('units','inch','position',[0,0,8,5]);
t=tiledlayout(2,3);



nexttile
colors=[0,0,0;0.5,0.5,0.5;0.231372549019608,0.298039215686275,0.752941176470588;0.705882352941177,0.0156862745098039,0.1490196078431370];
set(gca,'linestyleorder',{'-','-.','--',':'}, 'colororder',colors,'nextplot','add')

 set(0, 'DefaultFigureRenderer', 'painters');

Res=[500,400,300,200];
h=0.5;
l1=1.5;
l2=2;
lc=0.025;

for i=1:length(Res)
Re=Res(i);
path='data/Re_'+string(Re)+'h-'+string(h)+'l1-'+string(l1)+'l2-'+string(l2);
xbase=importdata(path+'/x_base.txt');
 sr=importdata(path+'/wss.txt');
 [xb,sr]=sort_v(xbase,sr);
  plot(xb,Re*sr); hold on

end
xlim([-7,7])
box on
ylabel('Scaled WSR, $Re \dot{\gamma}$','Interpreter','latex');
%  xlabel('$z$ ','Interpreter','latex');
title('  (a) Reynolds Number','Interpreter','latex')
legend('$Re=500$','$Re=400$','$Re=300$','$Re=200$')



nexttile
set(gca,'linestyleorder',{'-','-.','--',':'}, 'colororder',colors,'nextplot','add')
l2s=[2,3,4,5];
Re=500;
for i=1:length(l2s)
    l2=l2s(i);
path='data/Re_'+string(Re)+'h-'+string(h)+'l1-'+string(l1)+'l2-'+string(l2);
xbase=importdata(path+'/x_base.txt');
 sr=importdata(path+'/wss.txt');
 [xb,sr]=sort_v(xbase,sr);
  plot(xb,Re*sr); hold on
end
xlim([-7,7])
%  xlabel('$z$ ','Interpreter','latex');
legend('$\hat{l}_2=2$','$\hat{l}_2=3$','$\hat{l}_2=4$','$\hat{l}_2=5$')
box on
title('(b) Steepness','Interpreter','latex')

nexttile
set(gca,'linestyleorder',{'-','-.','--',':'}, 'colororder',colors,'nextplot','add')
hs=[0.5,0.4,0.3,0.2];
l2=2;
Re=500;
for i=1:length(hs)
    h=hs(i);
path='data/Re_'+string(Re)+'h-'+string(h)+'l1-'+string(l1)+'l2-'+string(l2);
xbase=importdata(path+'/x_base.txt');
 sr=importdata(path+'/wss.txt');
 [xb,sr]=sort_v(xbase,sr);
  plot(xb,Re*sr); hold on
end
xlim([-7,7])
legend('$\hat{h}=0.5$','$\hat{h}=0.5$','$\hat{h}=0.3$','$\hat{h}=0.2$')
box on
title('(c) Stenosis height','Interpreter','latex')


nexttile
l2=2;
h=0.5;

set(gca,'linestyleorder',{'-','-.','--',':'}, 'colororder',colors,'nextplot','add')
for i=1:length(Res)
Re=Res(i);
path='data/Re_'+string(Re)+'h-'+string(h)+'l1-'+string(l1)+'l2-'+string(l2);
xbase=importdata(path+'/x_base.txt');
 sr=importdata(path+'/base_ext.txt');
 [xb,sr]=sort_v(xbase,sr);
  plot(xb,sqrt(sr)); hold on
end
xlim([-7,7])
ylabel('VWF Extension, $\mathcal{E}$','Interpreter','latex');
ylim([0,22])
lk=yline(21.21);
lk.LineStyle='-.';
 xlabel('$z$ ','Interpreter','latex');
box on
nexttile
set(gca,'linestyleorder',{'-','-.','--',':'}, 'colororder',colors,'nextplot','add')

Re=500;
l2=2;
l1=1.5;
h=0.5;
for i=1:length(l2s)
    l2=l2s(i);
path='data/Re_'+string(Re)+'h-'+string(h)+'l1-'+string(l1)+'l2-'+string(l2);
xbase=importdata(path+'/x_base.txt');
 sr=importdata(path+'/base_ext.txt');
 [xb,sr]=sort_v(xbase,sr);
  plot(xb,sqrt(sr)); hold on

end
xlim([-7,7])
lk=yline(21.21);
lk.LineStyle='-.';
 xlabel('$z$ ','Interpreter','latex');
box on
ylim([0,22])

nexttile
set(gca,'linestyleorder',{'-','-.','--',':'}, 'colororder',colors,'nextplot','add')


hs=[0.5,0.4,0.3,0.2];
l2=2;
Re=500;
for i=1:length(hs)
    h=hs(i);
path='data/Re_'+string(Re)+'h-'+string(h)+'l1-'+string(l1)+'l2-'+string(l2);
xbase=importdata(path+'/x_base.txt');
 sr=importdata(path+'/base_ext.txt');
 [xb,sr]=sort_v(xbase,sr);
  plot(xb,sqrt(sr)); hold on
end
xlim([-7,7])
lk=yline(21.21);
lk.LineStyle='-.';
% legend('$h=0.5$','$h=0.5$','$h=0.3$','$h=0.2$')
box on
ylim([0,22])

t.TileSpacing = 'compact';
t.Padding = 'compact';
exportgraphics(figure1,'figs/total_dynamics.eps','ContentType','vector')


function[xb,var_sort]=sort_v(x,var)


A=[x,var];
B=sortrows(A,1);
xb=B(:,1);
var_sort=B(:,2:end);
end