%% Program for the 2D transport of a passive scalar
% Mateus Dias Ribeiro
% 24.05.14

clear all; % Clear workspace
clc; % Clear screen

% Defining the mesh density
%Ima = 64; Jma = 40;
Ima = 128; Jma = 80;
%Ima = 192; Jma = 120;
%Ima = 512; Jma = 320;

% Defintion of the grid variables
Ifim = 1; Ifi = 2; Ila = Ima+1; Ilap = Ima + 2;
Jfim = 1; Jfi = 2; Jla = Jma+1; Jlap = Jma + 2;

% l1 and l6 are the limits of the jet opening size
l1 = round(Jma/2-(0.05*Jma));
l2 = l1 + 1;
l3 = l2 + 1;
l6 = round(Jma/2+(0.05*Jma));
l5 = l6 - 1;
l4 = l5 -1;

%% Definition of variables
Lx = 0.5*Ilap/1000; % [m] Size of the domain x
Ly = 0.5*Jlap/1000; % [m] Size of the domain y
dx = Lx/Ilap; dy = Ly/Jlap; % [m] Size of the cell in x and y direction
V = dx*dy*dx; % [m^3] Volume of the cell in cubic meters
ni = 0.65*10^-5; % [m^2/s] Kinematic viscosity
Dx = ni; Dy = ni; % [m^2/s] Diffusion constant gamma in x and y direction
rho = 1.2; % [kg/m^3] Density
U = 0.8; % [m/s] Velocity of the jet
Re = U*Ly*((l6-l1)/Jlap)/ni; % Reynolds number
Ax = dy*dx; Ay = dx*dx; % [m^2] Area of the cell in x and y direction
CFL = 0.2; % CFL condition
dt = CFL*0.5*(dx+dy)/U; % [s] Time-step

%% Initialization of fields

% Initializing boundaries for u and v
u(1:Ilap,1:Jlap) = 0.0;
u(Ifim:Ilap,l1:l2) = 0.5*U;
u(Ifim:Ilap,l3:l4) = 1.0*U;
u(Ifim:Ilap,l5:l6) = 0.5*U;
v(1:Ilap,1:Jlap) =0.0;

phi(1:Ilap,1:Jlap) = 0.0;
phi(Ifim:Ilap,l1:l2) = 0.5*U;
phi(Ifim:Ilap,l3:l4) = 1.0*U;
phi(Ifim:Ilap,l5:l6) = 0.5*U;

% Initializing momentum in x and y direction
rhou = rho.*u;
rhov = rho.*v;
rhophi = rho.*phi;
[ue, vn] = mom2vel(rhou,rhov,rho);

%% Time-loop
cont1 = 0;
n = 0;
Nt = 120000;
for (cont1=1:Nt)
    dt = (CFL*0.5*(dx+dy))/(max(max(max(ue)),max(max(vn))));  % Dynamic Time-step
    Pe_max = (max(max(max(ue)),max(max(vn))))/(Dx/dx); % Max Peclet number
    %% Momentum    
    % Calculating convective fluxes in x and y directions for rhou
    [fluxConXu, fluxConYu] = calcFluxConCDS(rhou, Ax, Ay, ue, vn);
    
    % Calculating convective fluxes in x and y directions for rhov
    [fluxConXv, fluxConYv] = calcFluxConCDS(rhov, Ax, Ay, ue, vn);
    % Calculating diffusive fluxes in x and y directions for rhou
    [fluxDifXu, fluxDifYu] = calcFluxDifCDS(Ax, Ay, rhou, dx, dy, Dx, Dy);
    
    % Calculating diffusive fluxes in x and y directions for rhov
    [fluxDifXv, fluxDifYv] = calcFluxDifCDS(Ax, Ay, rhov, dx, dy, Dx, Dy);
    % Apply the fluxes for rhou and rhov
    rhou(Ifi:Ila,Jfi:Jla) = rhou(Ifi:Ila,Jfi:Jla) + dt./(V).*(fluxConXu(Ifim:Ila-1,Jfi:Jla)-fluxConXu(Ifi:Ila,Jfi:Jla)...
                                                            +fluxConYu(Ifi:Ila,Jfim:Jla-1)-fluxConYu(Ifi:Ila,Jfi:Jla)...
                                                            +fluxDifXu(Ifi:Ila,Jfi:Jla)+fluxDifYu(Ifi:Ila,Jfi:Jla));
                                              
    rhov(Ifi:Ila,Jfi:Jla) = rhov(Ifi:Ila,Jfi:Jla) + dt./(V).*(fluxConXv(Ifim:Ila-1,Jfi:Jla)-fluxConXv(Ifi:Ila,Jfi:Jla)...
                                                            +fluxConYv(Ifi:Ila,Jfim:Jla-1)-fluxConYv(Ifi:Ila,Jfi:Jla)...
                                                            +fluxDifXv(Ifi:Ila,Jfi:Jla)+fluxDifYv(Ifi:Ila,Jfi:Jla));  
   
                                              
    %% Scalar transport
    % Calculating convective fluxes in x and y directions for rhophi
    [fluxConXphi, fluxConYphi] = calcFluxConCDS(rhophi, Ax, Ay, ue, vn);
    
    % Calculating diffusive fluxes in x and y directions for rhophi
    [fluxDifXphi, fluxDifYphi] = calcFluxDifCDS(Ax, Ay, rhophi, dx, dy, Dx, Dy);
    
    % Apply the fluxes for rhophi
    rhophi(Ifi:Ila,Jfi:Jla) = rhophi(Ifi:Ila,Jfi:Jla) + dt./(V).*(fluxConXphi(Ifim:Ila-1,Jfi:Jla)-fluxConXphi(Ifi:Ila,Jfi:Jla)...
                                                            +fluxConYphi(Ifi:Ila,Jfim:Jla-1)-fluxConYphi(Ifi:Ila,Jfi:Jla)...
                                                            +fluxDifXphi(Ifi:Ila,Jfi:Jla)+fluxDifYphi(Ifi:Ila,Jfi:Jla));
  
    
    %% Boundary treatment
    rhou(Ilap,:) = max(0,rhou(Ila,:)); %inhibits Inflow
    rhou(:,Jfim) = rhou(:,Jfi);          
    rhou(:,Jlap) = rhou(:,Jla);
     
    rhov(Ilap,:) = max(0,rhov(Ila,:)); %inhibits Inflow
    rhov(Ifim,:) = U.*(rand(1,Jlap)-0.5); % Generates turbulence at inlet
    %rhov(Ifim,:) = rhov(Ifi,:);
    rhov(:,Jfim) = rhov(:,Jfi);   
     
    rhophi(Ilap,:) = rhophi(Ila,:); %inhibits Inflow
    rhophi(:,Jfim) = rhophi(:,Jfi);          
    rhophi(:,Jlap) = rhophi(:,Jla);
    
    [ue, vn] = mom2vel(rhou,rhov,rho);
    
    %% Poisson equation
    al=1;
    % Calculate divergence
    [div]=calcDivergence(ue,vn,dx,Ilap,Jlap);
    % Calculate the potential field PHI
    [PHI]=PoissonSolver(div,Ilap,Jlap,dx);
    
    %% Pressure correction for rhou and rhov
    ue(1:Ilap-1,1:Jlap)=ue(1:Ilap-1,1:Jlap) + al*(PHI(2:Ilap,1:Jlap)-PHI(1:Ilap-1,1:Jlap));
    vn(1:Ilap,1:Jlap-1)=vn(1:Ilap,1:Jlap-1) + al*(PHI(1:Ilap,2:Jlap)-PHI(1:Ilap,1:Jlap-1));
    
    rhou(2:Ilap-1,2:Jlap-1)=rhou(2:Ilap-1,2:Jlap-1) + al*rho/2*(PHI(3:Ilap,2:Jlap-1)-PHI(1:Ilap-2,2:Jlap-1));
    rhov(2:Ilap-1,2:Jlap-1)=rhov(2:Ilap-1,2:Jlap-1) + al*rho/2*(PHI(2:Ilap-1,3:Jlap)-PHI(2:Ilap-1,1:Jlap-2));
        
    u = rhou./rho;
    v = rhov./rho;
    phi = rhophi./rho;
    
    %% Plot
    if (mod(cont1,100) == 0); pause(1.0); 
        imagesc(phi');
        %imagesc(PHI');
        %contourf(phi'), grid on;
        %imagesc(PHI);
    title({['Scalar transport, Re = ',num2str(Re)];
            ['Time-step = ',num2str(cont1)]
            ['dt = ',num2str(dt)]
            ['Pe = ',num2str(Pe_max)]})
    end
    n = n + dt;
    cont1 = cont1 + 1;
end