function [fluxDifX, fluxDifY] = calcFluxDifCDS(Ax, Ay, rhophi, dx, dy, Dx, Dy)

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

fluxDifX = zeros(I,J);
fluxDifY = zeros(I,J);

fluxDifX(Ifi:Ila,:) = (Dx*Ax/dx).*(rhophi(Ifim:Ila-1,:)-2.*rhophi(Ifi:Ila,:)+rhophi(Ifi+1:Ilap,:));
fluxDifY(:,Jfi:Jla) = (Dy*Ay/dy).*(rhophi(:,Jfim:Jla-1)-2.*rhophi(:,Jfi:Jla)+rhophi(:,Jfi+1:Jlap));