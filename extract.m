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