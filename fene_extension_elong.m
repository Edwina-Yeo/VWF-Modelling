

function [length_v,A]=fene_extension_elong(params,if_mod,sr_vec)

hatb=params(2);
gamma_star=params(3);
De=params(1);
LL=params(end);
delta=params(4);
N=length(sr_vec);

A=zeros(N,2);

length_v=zeros(N,1);
previous=[1,1];
for i=1:length(sr_vec)

DATA = [hatb,gamma_star,De,LL,delta,sr_vec(i)];
options = optimset('Display','off');% so that fsolve doesnt print

% f = @(y) fene_shear(y,DATA); 

f = @(y) fene_elong(y,DATA,if_mod); 
sol = fsolve(f, previous,options);
previous=sol;
length_v(i)=sqrt((sol(1)+sol(2)+1)/3);
A(i,1)=sol(1);
A(i,2)=sol(2);
end
