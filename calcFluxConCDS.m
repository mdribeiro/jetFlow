function [fluxConX, fluxConY] = calcFluxConCDS(rhophi, Ax, Ay, ue, vn)

[I, J] = size(rhophi); 
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

fluxConX=zeros(I,J);
fluxConY=zeros(I,J);

fluxConX(Ifim:Ila,:) = (Ax/4).*ue(Ifim:Ila,:).*(rhophi(Ifim:Ila,:)+rhophi(Ifi:Ilap,:));
fluxConY(:,Jfim:Jla) = (Ay/4).*vn(:,Jfim:Jla).*(rhophi(:,Jfim:Jla)+rhophi(:,Jfi:Jlap));
