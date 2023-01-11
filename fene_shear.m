
% Algebraic system for the extension of VWF in 2D shear flow
function [F,J] = fene_shear(A,data)
% A=[Axx,Axy,Ayy]
hatb=data(1);
gamma_star=data(2);
De=data(3);
LL=data(4);
delta=data(5);
sr=data(6);
tau=De*((tanh(hatb*(sr-gamma_star))+1)/2+delta);

f=LL^2/(LL^2-A(1)-A(3));
a=LL^2/(LL^2-2);
F(1) = 2*(sr*tau)*A(2)-(A(1)*f-a);
F(2)=(tau*sr)*A(3)-(A(2)*f);
F(3)=(A(3)*f-a);
% if nargout > 1
% % form jacobian
% J=zeros(3,3);
% J(1,1)=-f-LL^2*A(1)/(LL^2-A(1)-A(3)-1)^2;
% J(1,2)=2*(sr*tau);
% J(1,3)=-LL^2*A(1)/(LL^2-A(1)-A(3)-1)^2;
% 
% J(2,1)=-LL^2*A(2)/(LL^2-A(1)-A(3)-1)^2;
% 
% J(2,2)=-f;
% 
% J(2,3)=(tau*sr)-LL^2*A(2)/(LL^2-A(1)-A(3)-1)^2;
% 
% 
% J(3,1)=-LL^2*A(3)/(LL^2-A(1)-A(3)-1)^2;
% J(3,3)=f+LL^2*A(3)/(LL^2-A(1)-A(3)-1)^2;
% 
% end
end

