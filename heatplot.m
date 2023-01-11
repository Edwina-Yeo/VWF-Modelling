%%



% step 1: make a continuous line for the recirculation zone
[X,Y]=meshgrid(Res,xb); %for plotting
z_list=linspace(-l1-l2,l1+l2+0.2,200);
re_list=linspace(200,500,200);
[Re_i,z_i]=meshgrid(re_list,z_list);
BDmatrixq = interp2(X,Y,SS,Re_i,z_i,'spline');
cont=contour(z_i,Re_i,BDmatrixq,[0,100],'k','LineWidth',1.5); hold on

BDmatrixq = interp2(X,Y,SS,Re_i,z_i,'spline');

% load colormap
load('red_blue_cmap.mat')
tdata = load('red_blue_cmap.mat');  %load into structure

figure1=figure('units','inch','position',[0,0,6,3]);

t=tiledlayout(1,2)
nexttile
heat_plot_Re(abs(WSS),'(a) Wall Shear Rate',cont,tdata,X,Y,l1,l2,path_outer)
nexttile

% heat_plot_Re(tau/0.04,'(b) VWF Relaxation Time',cont,tdata,X,Y,l1,l2,path_outer)

% colormap(tdata.cmap)  
plot3(cont(1,2:end),cont(2,2:end),10*ones(length(cont(2,2:end)),1),'k--','LineWidth',1.5); hold on
p=xline(0,'k');
p.LineStyle='-';
p.LineWidth=1.5;


view(2)
xlim([-l1-l2,l1+l2+0.2])
colorbar('TickLabelInterpreter','latex');
ylim([200,500])
ylabel('$Re$','Interpreter','latex');
xlabel('$z$ ','Interpreter','latex');
title('(b) VWF Relaxation Time','Interpreter','latex')

heat_plot_Re(sqrt(EXT+0.01),'(b) VWF Extension',cont,tdata,X,Y,l1,l2,path_outer)

t.Padding='compact'
exportgraphics(figure1,'meta.png','Resolution',300) 


 
function c = redblue(m)
%REDBLUE    Shades of red and blue color map
%   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with bright blue, range through shades of
%   blue to white, and then through shades of red to bright red.
%   REDBLUE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(redblue)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.
%   Adam Auton, 9th October 2009
if nargin < 1, m = size(get(gcf,'colormap'),1); end
if (mod(m,2) == 0)
    % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
    m1 = m*0.5;
    r = (0:m1-1)'/max(m1-1,1);
    g = r;
    r = [r; ones(m1,1)];
    g = [g; flipud(g)];
    b = flipud(r);
else
    % From [0 0 1] to [1 1 1] to [1 0 0];
    m1 = floor(m*0.5);
    r = (0:m1-1)'/max(m1,1);
    g = r;
    r = [r; ones(m1+1,1)];
    g = [g; 1; flipud(g)];
    b = flipud(r);
end
c = [r g b]; 

end



function cmp = getPyPlot_cMap(nam,n,keepAlpha,pyCmd)
% cmp = getPyPlot_cMap(nam [, n, keepAlpha, pyCmd])
%
%
% ::INPUT::
%
% nam:      Colormap name from matplotlib library (case sensitve!). See
%           below. Alternatively, '!GetNames' returns a cellstring of
%           available colormaps.
% n:        Number of Colors; defaults to 128
% keepAlpha: Switch to keep the alpha channel of the colormap (4th colum);
%           defaults to false. If true, a Nx4 matrix is returned in cmp
%           (instead of Nx3).
% pyCmd:    python command; defaults to 'python'
%
% 
% ::OUTPUT::
%
% cmp       A Nx3 (Nx4 if keepAlpha is true) double array of RGB(A) values.
%
%
% Colormap name can be one of the following:
%
%     Accent      gist_earth        Oranges           RdYlBu
%     Accent_r    gist_earth_r      Oranges_r         RdYlBu_r
%     afmhot      gist_gray         OrRd              RdYlGn
%     afmhot_r    gist_gray_r       OrRd_r            RdYlGn_r
%     autumn      gist_heat         Paired            Reds
%     autumn_r    gist_heat_r       Paired_r          Reds_r
%     binary      gist_ncar         Pastel1           seismic
%     binary_r    gist_ncar_r       Pastel1_r         seismic_r
%     Blues       gist_rainbow      Pastel2           Set1
%     Blues_r     gist_rainbow_r    Pastel2_r         Set1_r
%     bone        gist_stern        pink              Set2
%     bone_r      gist_stern_r      pink_r            Set2_r
%     BrBG        gist_yarg         PiYG              Set3
%     BrBG_r      gist_yarg_r       PiYG_r            Set3_r
%     brg         GnBu              PRGn              Spectral
%     brg_r       GnBu_r            PRGn_r            spectral
%     BuGn        gnuplot           prism             spectral_r
%     BuGn_r      gnuplot_r         prism_r           Spectral_r
%     BuPu        gnuplot2          PuBu              spring
%     BuPu_r      gnuplot2_r        PuBu_r            spring_r
%     bwr         gray              PuBuGn            summer
%     bwr_r       gray_r            PuBuGn_r          summer_r
%     CMRmap      Greens            PuOr              terrain
%     CMRmap_r    Greens_r          PuOr_r            terrain_r
%     cool        Greys             PuRd              winter
%     cool_r      Greys_r           PuRd_r            winter_r
%     coolwarm    hot               Purples           Wistia
%     coolwarm_r  hot_r             Purples_r         Wistia_r
%     copper      hsv               rainbow           YlGn
%     copper_r    hsv_r             rainbow_r         YlGn_r
%     cubehelix   jet               RdBu              YlGnBu
%     cubehelix_r jet_r             RdBu_r            YlGnBu_r
%     Dark2       nipy_spectral     RdGy              YlOrBr
%     Dark2_r     nipy_spectral_r   RdGy_r            YlOrBr_r
%     flag        ocean             RdPu              YlOrRd
%     flag_r      ocean_r           RdPu_r            YlOrRd_r
% 
% V 1.3; Konrad Schumacher, 07.2019
if strcmpi(nam,'!GetNames')
    % translate switch to retrieve colormap names into python-switch:
    nam = 'listCMapNames';
end
% defaults:
if ~exist('n','var') || isempty(n)
    n = 128;
end
if ~exist('keepAlpha','var') || isempty(keepAlpha)
    keepAlpha = 0;
end
if ~exist('pyCmd','var') || isempty(pyCmd)
    pyCmd = 'python';
end
% check if python script is present
pyScript = which('pyplotCMap2txt.py');
assert(~isempty(pyScript), 'getPyPlot_cMap:PyScriptNotFound', ...
    'Could not find python script (%s).','pyplotCMap2txt.py');
tmpf = tempname;
% call python script
comd = sprintf('%s "%s" %s -o "%s" -n %d', pyCmd, pyScript, nam, tmpf, n);
[s,m] = system(comd);
% check if system command ran w/o error
assert(s==0, 'getPyPlot_cMap:SystemCMDFailed', ...
        'There was an error executing the command\n\t%s\nSystem returned:\n\t%s', ...
        comd, m);
if strcmp(nam,'listCMapNames')
    % cMap names retrieved; read text file
    fid = fopen(tmpf,'r');
    cmp = textscan(fid,'%s');
    fclose(fid);
    cmp = cmp{1};
    
else
    % load cMap data from text file
    cmp = load(tmpf,'-ascii');
    if keepAlpha
    else % remove 4th column of alpha values
        cmp = cmp(:,1:3);
    end
end
% delete temp-file
delete(tmpf);
end
function c = coolwarm(m)
%COOLWARM    Blue to white to red divergent color map
%   COOLWARM(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with dark blue, range through shades of
%   blue to white, and then through shades of red to dark red.
%   As seen in "Diverging Color Maps for Scientific Visualization"
%   by Kenneth Moreland of Sandia National Laboratories.
%
%   COOLWARM, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(coolwarm)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, COLORMAP, RGBPLOT.
%
%   Created by:   Sam Borgeson 
%   Email:        sborgeson@berkeley.edu 
%   Last update:  1/30/2012
%

if nargin < 1, m = size(get(gcf,'colormap'),1); end

%   Seed color values from "Diverging Color Maps for Scientific
%   Visualization" by Kenneth Moreland of Sandia National Laboratories
a = [   59 76 192;...
        68 90 204;...
        77 104 215;...
        87 117 225;...
        98 130 234;...
        108 142 241;...
        119 154 247;...
        130 165 251;...
        141 176 254;...
        152 185 255;...
        163 194 255;...
        174 201 253;...
        184 208 249;...
        194 213 244;...
        204 217 238;...
        213 219 230;...
        221 221 221;...
        229 216 209;...
        236 211 197;...
        241 204 185;...
        245 196 173;...
        247 187 160;...
        247 177 148;...
        247 166 135;...
        244 154 123;...
        241 141 111;...
        236 127 99;...
        229 112 88;...
        222 96 77;...
        213 80 66;...
        203 62 56;...
        192 40 47;...
        180 4 38] / 255;

% get the number of color entries (= rows of a)
n = size(a,1);
% interpolate the originals to have m values as requested
% note that there are n-1 steps between 1 and n
c = interp1(1:n,a,1:(n-1)/(m-1):n);
end

function heat_plot_Re(var,var_name,cont,tdata,X,Y,l1,l2,path_outer)
% figure1=figure('units','inch','position',[0,0,3,3]);

colormap(tdata.cmap)  
plot3(cont(1,2:end),cont(2,2:end),10*ones(length(cont(2,2:end)),1),'k--','LineWidth',1.5); hold on
p=xline(0,'k');
p.LineStyle='-';
p.LineWidth=1.5;

k=surf(Y,X,var);hold on

k.EdgeColor = 'none';
view(2)
xlim([-l1-l2,l1+l2+0.2])
colorbar('TickLabelInterpreter','latex');
ylim([200,500])
ylabel('$Re$','Interpreter','latex');
xlabel('$z$ ','Interpreter','latex');
title(var_name,'Interpreter','latex')
% exportgraphics(figure1,append(path_outer,'/',var_name,'.png'),'Resolution',300) 

end
