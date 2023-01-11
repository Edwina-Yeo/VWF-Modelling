%

figure
    [x1,y1] = meshgrid(-2:.2:2,-1:.15:1);
    z1 = x1 .* exp(-x1.^2 - y1.^2);
    [u,v,w] = surfnorm(x1,y1,z1);
    q = quiver3(x1,y1,z1,u,v,w); hold on; surf(x1,y1,z1); hold off;
    drawnow;                  % This is needed as, if the plot is not plotted, the VertexData for the quiver are not existent.
    mag = 0.1*ones(size(u));  % Length of the vectors : 0.1
    
    size(mag)
    SetQuiverLength(q,mag,'HeadLength',0.05,'HeadAngle',90);


%%

%     set(0, 'DefaultFigureRenderer', 'opengl');
% 
% figure1=figure('units','inch','position',[0,0,7,2]);
% ax_vel = axes;
% s=trisurf(DT(IO, :),x,y,-ztot,'FaceColor','interp'); hold on
% s.EdgeColor='none';
% view(2)
% xlim([-10,30]);
% plot3([-10,30,30,-10],[0,0,1,1],[100,100,100,100],'k','LineWidth',0.2)
% ylabel(ax_vel,'$r$','Interpreter','latex')
% a=colorbar(ax_vel,'TickLabelInterpreter','latex');
% a.Label.String = ' $\dot{\gamma}-\dot{\omega}$';
% a.Label.Interpreter = 'latex';
% % a.Limits=[-1.92,1.92]
% colormap(tdata.cmap)
% 
% title(ax_vel,'Flow classification','Interpreter','latex')


[xq,yq] = meshgrid(-10:0.5:30, 0.05:0.1:0.95);
vr = griddata(x,y,zr,xq,yq);
vz = griddata(x,y,zz,xq,yq);
polyin = polyshape(x(k),y(k));

xflat=reshape(xq,1,[]);
yflat=reshape(yq,1,[]);
[TFin,TFon] = isinterior(polyin,xflat,yflat);
indIN=reshape(TFin,size(xq));
indon=reshape(TFon,size(xq));
indon=not(indon);





AR=vr.*indIN.*indon;
AR(AR==0)=nan;
% % AR(AR==nan)=0;
AZ=vz.*indIN.*indon;
AZ(AZ==0)=nan;
% AZ(AZ==nan)=0;


rr=sqrt(abs(AR)).*sign(AR)-1
rz=sqrt(abs(AZ)).*sign(AZ)-1
Ln = sqrt(rr.^2+rz.^2)+1;


log(rr)
figure
quiver(xq,yq,rz./Ln.*log(Ln),-rr./Ln.*log(Ln)); hold on


%%

% figure1=figure('units','inch','position',[0,0,7,2]);
% ax_vel = axes;
% s=trisurf(DT(IO, :),x,y,log(zsr+0.26),'FaceColor','interp'); hold on
% s.EdgeColor='none';
% view(2)
% xlim([-10,30]);
% plot3([-10,30,30,-10],[0,0,1,1],[100,100,100,100],'k','LineWidth',0.2)
% ylabel(ax_vel,'$r$','Interpreter','latex')
% a=colorbar(ax_vel,'TickLabelInterpreter','latex');
% a.Label.String = ' $\dot{\gamma}-\dot{\omega}$';
% a.Label.Interpreter = 'latex';
% % a.Limits=[-1.92,1.92]
% colormap(tdata.cmap)
% 
% title(ax_vel,'Flow classification','Interpreter','latex')


ztot=abs(zsr)-abs(zrot)

figure1=figure('units','inch','position',[0,0,7,2]);
ax_vel = axes;
s=trisurf(DT(IO, :),x,y,ztot,'FaceColor','interp'); hold on
s.EdgeColor='none';
view(2)
xlim([-10,30]);
plot3([-10,30,30,-10],[0,0,1,1],[100,100,100,100],'k','LineWidth',0.2)
ylabel(ax_vel,'$r$','Interpreter','latex')
a=colorbar(ax_vel,'TickLabelInterpreter','latex');
a.Label.String = ' $\dot{\gamma}-\dot{\omega}$';
a.Label.Interpreter = 'latex';
a.CAxis=[-1.92,1.92]
colormap(tdata.cmap)

title(ax_vel,'Flow classification','Interpreter','latex')

% 
% 
% figure1=figure('units','inch','position',[0,0,7,2]);
% ax_vel = axes;
% s=trisurf(DT(IO, :),x,y,log(sqrt(ze+0.05)),'FaceColor','interp'); hold on
% s.EdgeColor='none';
% view(2)
% xlim([-10,30]);
% plot3([-10,30,30,-10],[0,0,1,1],[100,100,100,100],'k','LineWidth',0.2)
% ylabel(ax_vel,'$r$','Interpreter','latex')
% a=colorbar(ax_vel,'TickLabelInterpreter','latex');
% a.Label.String = ' $\dot{\gamma}-\dot{\omega}$';
% a.Label.Interpreter = 'latex';
% % a.Limits=[-1.92,1.92]
% colormap(tdata.cmap)
% 
% title(ax_vel,'Flow classification','Interpreter','latex')

% qH2=quiver3(xq,yq,2*ones(size(xq)),AZ_norm,AR_norm,W); hold on
 qH2.ShowArrowHead = 'off';
 qH2.Color = 'k';
% % %     drawnow; 
quiver3(xq,yq,ones(size(xq)),rz./Ln.*log(Ln),-rr./Ln.*log(Ln),zeros(size(xq))); hold on

%     mags = [log(sqrt(abs(AZ))+1),0*log(sqrt(abs(AR))+1),zeros*ones(size(xq))];  % Length of the vectors : 0.1
%  SetQuiverLength(qH2,mags,'HeadLength',0.1,'HeadAngle',90); 

%  ylim([0,1])
% xlabel('$z$','Interpreter','latex')
% ylabel('$r$','Interpreter','latex') 
% grid off
% view(2)
% 
% colormap(tdata.cmap) 
% % exportgraphics(figure1,'data/orient4.png','Resolution',300) 
% plot3([-10,30,30,-10,-10],[0,0,1,1,0],[100,100,100,100,100],'k','LineWidth',0.2)
% box on
% exportgraphics(figure1,'data/orient2.png','Resolution',300) 
% function SetQuiverLength(q,mags,varargin)
%--------------------------------------------------
% function SetQuiverLength(q,mags)
%
%   Function that sets the length of the quiver
%   You can give a list of length (in x,y axis units), the function will rescale the vectors to the right length
%
% Input: q    = handle to quiver plot, can be quiver or quiver3
%        mags = Desired length of each vector in units of the x,y axis
%               If mags is a scalar, all the vector will have that length
% Optional inputs: given as ('variable',value)
%        'HeadLength' = single value for the length of the head in units of the xyz axis.
%                       It is applied to all the vectors
%        'HeadAngle' = angle between the two lines forming the head
%                      default=28.0724^\circ
%        'RotHead'   = Angle [deg] by which the head will be rotated around the vector axis.
%                      This allows to set the head in different planes
%
% NOTE: For some unknown reason, MATLAB does not always simply change the length of the vectors and requires a
%       pause statement towards the end.
%       Would your vectors not be the right size, increase the duration of the pause
%       (or even better, suggest a solution that does not require a pause)
%
%%
close all
figure
% subplot(1,2,1)
% scatter(ztot,ze,[],zsr); hold on
% subplot(1,2,2)
colorbar
zrot=ztot+zsr;
ztot=abs(zsr)-abs(zrot)

plot(x,ztot,'.'); hold on

plot(x,zrot,'.'); hold on
plot(x,zsr,'.')
%%
figure
scatter3(ztot,log(ze),x,[],log(ze)); hold on
view(3)
colorbar

%%

figure
scatter(ztot,zsr)

[xq,yq] = meshgrid(min(ztot):.2:max(ztot), min(zsr):1:max(zsr));
vq = griddata(ztot,zsr,ze,xq,yq);

%%

xdata=ztot.*abs(zsr)

data = [xdata,log(ze-min(ze)+0.01)];


figure
[hh3,N] = hist3(data, 'Nbins',[1,1]*200);
% figure
image(flipud(hh3))

% surf(% ax = gc
[X,Y]=meshgrid(N{1},N{2})
s=surf(X,Y,log(hh3+0.01)')
view(2)
s.EdgeColor='none'
% xt = ax.XTick;
% xt=ztot;
% % yt = ax.YTick;
% yt=zsr;
% view(3)
% ax.XTickLabel = xt*10;
% set(ax, 'YTick',[0 yt], 'YTickLabel', [flip([0 yt])]*10)

% figure
% sols=[zsr,ztot,ze];
% 
% plot3(zsr,ztot,ze,'.')
% xlabel('sr')
% ylabel('tot')
% zlabel('ext')


% biplot(coeff(:,1:2),'scores',score(:,1:2),'varlabels',{'v_1','v_2','v_3','v_4'});

% plot(zsr,ztot,'.'); hold on
