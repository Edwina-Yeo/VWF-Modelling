
function F = fene_elong(A,data,if_mod)


hatb=data(1);
gamma_star=data(2);
De=data(3);
L=data(4);
delta=data(5);
sr=data(6);
if if_mod==1
tau=De*((tanh(hatb*(sr-gamma_star))+1)/2+delta);
else
    tau=De;
end
% F(1) = (sr*tau)*A(1)*(sqrt(2))-L^2*(A(1)/(L^2-A(1)-A(2))-1/(L^2-3));
% F(2)=(sr*tau)*A(2)*(sqrt(2))+L^2*(A(2)/(L^2-A(1)-A(2))-1/(L^2-3));
F(1) = (sr*tau)*A(1)*(2)-L^2*(A(1)/(L^2-A(1)-A(2))-1/(L^2-2));
F(2)=(sr*tau)*A(2)*(2)+L^2*(A(2)/(L^2-A(1)-A(2))-1/(L^2-2));
end