function [ue, vn] = mom2vel(rhou,rhov,rho)

[I, J] = size(rhou); 
Ima = I-2;
Jma = J-2;
Ifim = 1;
Jfim = 1;
Ifi = 2;
Ila = Ima+1;
Jfi = 2;
Jla = Jma+1;
Ilap = Ima+2;
Jlap = Jma+2;

ue(Ifim:Ila,Jfim:Jlap) = (rhou(Ifim:Ila,Jfim:Jlap)+rhou(Ifi:Ilap,Jfim:Jlap))/(rho*2);
vn(Ifim:Ilap,Jfim:Jla) = (rhov(Ifim:Ilap,Jfim:Jla)+rhov(Ifim:Ilap,Jfi:Jlap))/(rho*2);
