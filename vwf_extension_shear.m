
function [length_v,A]=vwf_extension_shear(params,sr_vals)
hatb=params(2);
gamma_star=params(3);
De=params(1);
LL=params(end);
delta=params(4);

N=length(sr_vals);
length_v=zeros(N,1);
previous=[1,0,1];
A=zeros(N,3);
for i=1:length(sr_vals)
DATA = [hatb,gamma_star,De,LL,delta,sr_vals(i)];
options = optimset('Display','off');% so that fsolve doesnt print
f = @(y) fene_shear(y,DATA); 
sol = fsolve(f, previous,options);
previous=sol;
length_v(i)=sqrt((sol(1)+sol(3)+1)/3);
A(i,1)=sol(1);
A(i,2)=sol(2);
A(i,3)=sol(3);
end
end 