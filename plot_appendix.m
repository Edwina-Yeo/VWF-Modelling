%% This generates the convergence plot for the the appendix
 

%---------------Import data from python----------------------------

mine=load('data/Convergence/min_extension.txt');
maxe=load('data/Convergence/max_extension.txt');
maxw=load('data/Convergence/max_axial_vel.txt');
verts=load('data/Convergence/mesh_vertices.txt');

%--------- plotting macros------------------------
set(groot,'DefaultAxesFontSize',10);
set(groot,'DefaultTextInterpreter','latex')
set(groot,'DefaultLegendInterpreter','latex')
set(groot,'DefaultAxesTickLabelInterpreter','latex')
set(groot,'DefaultTextFontname', 'CMU Serif')
set(groot,'DefaultAxesFontName', 'CMU Serif')
set(groot,'DefaultLineLineWidth',1)


%---------------Plot----------------------------

loc=4;% location of selected mesh 


figure1=figure('units','inch','position',[0,0,6,2]);
t=tiledlayout(1,1);
nexttile
 loglog(verts(2:end),abs((maxe(2:end)-maxe(1)))/(maxe(1)),'k.-','MarkerSize',10) ;hold on
loglog(verts(2:end),abs((maxw(2:end)-maxw(1)))/(maxw(1)),'.-','MarkerSize',10,'Color',[0.5,0.5,0.5]) ;hold on
 loglog(verts(loc),abs((maxe(loc)-maxe(1)))/(maxe(1)),'r.-','MarkerSize',10) ;hold on
loglog(verts(loc),abs((maxw(loc)-maxw(1)))/(maxw(1)),'r.-','MarkerSize',10) ;hold on


xlabel('Mesh vertices')
ylabel('Relative error')

xticks([1e2,1e3,1e4,1e5])
yticks([1e-5,1e-4,1e-3,1e-2,1e-1])

legend('VWF Eq.','Navier Stokes Eq.','Selected mesh','Location','southwest')

t.TileSpacing = 'compact';
t.Padding = 'compact';
exportgraphics(figure1,'data/convergence.eps','ContentType','vector')
