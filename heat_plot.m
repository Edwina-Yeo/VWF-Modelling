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

