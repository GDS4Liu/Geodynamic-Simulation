% fully staggered grid to solve 2D stokes equation
% using the Vx, Vz, P formulation (constant viscosity)
% Based on Boris Lecture
% 2023.04.27

%%
clear;

tic;

% mu = 1;
% Material properties
%     phase1 phase2
nu  = [1e20      1e22];  %Pa s
rho = [3200      3300];  %kg/m^3
Cp  = [1000      1100];  %J/kg/K
k   = [10        3   ];  %W/m/K
Ttop= 273; %K
Tbot= 1673;%K
T0  = 1600;%K
R   = 8.314; %J/(mol*K)
% k   = 2.2; % Come from mantle convection
Ea  = 5.4e+5;
Va  = 1.5e-5;
timesum  = 1e+5;

% dislocation creep
n   = 2.3;
F2  = 1/(2^((n-1)/n)*3^((n+1)/2*n));
m   = 0;
h   = 3e-3;
%AD  = 1./nu;

xnum = 61;
znum = 51;

xsize = 1000e+3;
zsize = 1000e+3;
g     = 10;

x_int = 0:xsize/100:xsize;
z_int = cos(x_int*2*pi/xsize)*zsize/10 -zsize/2 ;


% Setup the model
% Block
% xwidth = [0.4,0.6];
% zwidth = [-0.6,-0.4];


% model grid
dx = xsize/(xnum-1);
dz = zsize/(znum-1);
[xgrid,zgrid] = meshgrid(0:dx:xsize,-zsize:dz:0);


% Vx-grid
XVx = xgrid(1:end-1,:);
ZVx = (zgrid(1:end-1,:)+zgrid(2:end,:))*0.5;

% Vz-grid 
XVz = (xgrid(:,1:end-1)+xgrid(:,2:end))*0.5;
ZVz = zgrid(:,1:end-1);

% P-grid
XP  = (xgrid(1:end-1,1:end-1)+xgrid(2:end,2:end))*0.5;
ZP  = (zgrid(1:end-1,1:end-1)+zgrid(2:end,2:end))*0.5;


%% set the temperature grid
% Temperature in shear viscosity(on node)
Ts  = linspace(Tbot,Ttop,znum)';
% Temperature in centres of cell
Tn  = interp1(zgrid(:,1),Ts,ZP(:,1));
Ts  = repmat(Ts,1,xnum);
Tn  = repmat(Tn,1,xnum-1);

Rho = ones(znum,xnum)*rho(2);   % used for temperature equation
Rhos= ones(znum,xnum-1)*rho(2); % used for stokes equation
mu_n= ones(znum-1,xnum-1)*nu(2).*exp(1./(R*Tn)-1/(R*T0));
mu_s= ones(znum,xnum)*nu(2).*exp(1./(R*Ts)-1/(R*T0));
k_n = ones(znum,xnum)*k(2);
Cp_n= ones(znum,xnum)*Cp(2);

z_int_intp = interp1(x_int,z_int,XP(1,:));


for ix = 1:length(z_int_intp)
    ind = find(ZVz(:,1)<z_int_intp(ix));
    Rho(ind(1:end-1),ix)  = rho(1);
    Rhos(ind(1:end-1),ix) = rho(1);
    mu_n(ind(1:end-1),ix) = nu(1).*exp(1./(R*Tn(ind(1:end-1),ix))-1/(R*T0));
    mu_s(ind(1:end-1),ix) = nu(1).*exp(1./(R*Ts(ind(1:end-1),ix))-1/(R*T0));
    k_n(ind(1:end-1),ix)  = k(1);
    Cp_n(ind(1:end-1),ix) = Cp_n(1);
    
    
    fac = (z_int_intp(ix) - ZVz(ind(end),1))/dz;
    Rho(ind(1:end),ix)  = fac*rho(1) + (1-fac)*rho(2);
    Rhos(ind(end),ix)   = fac*rho(1) + (1-fac)*rho(2);
    mu_n(ind(end),ix)  = fac*nu( 1).*exp(1./(R*Tn(ind(end),ix))-1/(R*T0)) + (1-fac)*nu( 2).*exp(1./(R*Tn(ind(end),ix))-1/(R*T0));
    mu_s(ind(end),ix)  = fac*nu( 1).*exp(1./(R*Ts(ind(end),ix))-1/(R*T0)) + (1-fac)*nu( 2).*exp(1./(R*Ts(ind(end),ix))-1/(R*T0));
    k_n(ind(end),ix)   = fac*k ( 1) + (1-fac)*k(2);
    Cp_n(ind(end),ix)  = fac*Cp_n(1)+ (1-fac)*Cp_n(2);
    
end

% for ix = 1:length(z_int_intp)
%     ind = find(ZVz(:,1)<z_int_intp(ix));
%     Nu(ind(1:end-1),ix) = nu(1);
%     Rho(ind(1:end-1),ix) = rho(1);
%     
%     fac = (z_int_intp(ix) - ZVz(ind(end),1))/dz;
%     Rho(ind(end),ix) = fac*rho(1) + (1-fac)*rho(2);
%     Nu(ind(end),ix)  = fac*nu( 2) + (1-fac)*nu( 2);
%     
% end

% 
% for ix = 1:length(Px)
%     indz = find(ZVz(:,1)>zwidth(1) & ZVz(:,1)<zwidth(2) );
%     indx = find(XVz(1,:)>xwidth(1) & XVz(1,:)<xwidth(2));
%     Nu(indz(:),indx(:)) = nu(1).*exp(1./(R*T(indz(:),indx(:)))-1/(R*T0));
%     Rho(indz(:),indx(:)) = rho(1);
% %     
% %     fac = (ZVz(zwidth(2),:) - ZVz(indz(end),1))/dz;
% %     Rho(indz(end),ix) = fac*rho(1) + (1-fac)*rho(2);
% %     Mu(indz(end),ix)  = fac*nu( 2) + (1-fac)*nu( 2);
% 
% 
% end



Number_Phase = zeros(znum    +    znum -1, xnum + xnum-1);
Number_ind   = zeros(znum    +    znum -1, xnum + xnum-1);
Number_Vx    = zeros(znum-1 ,xnum  );
Number_Vz    = zeros(znum   ,xnum-1);
Number_P     = zeros(znum-1 ,xnum-1);
Number_T     = zeros(znum   ,xnum  );
for ix = 1:2:xnum+xnum-1, for iz = 2:2:znum+znum-1 , Number_Phase(iz,ix) = 1; end; end
for ix = 2:2:xnum+xnum-1, for iz = 1:2:znum+znum-1 , Number_Phase(iz,ix) = 2; end; end
for ix = 2:2:xnum+xnum-1, for iz = 2:2:znum+znum-1 , Number_Phase(iz,ix) = 3; end; end
num    = 1;
for ix = 1:size(Number_Phase,2)
    for iz = 1:size(Number_Phase,1)
        if Number_Phase(iz,ix) == 1
            Number_ind(iz,ix) = num;
            num = num+1;
        end
    end
end


for ix = 1:size(Number_Phase,2)
    for iz = 1:size(Number_Phase,1)
        if Number_Phase(iz,ix) == 2
            Number_ind(iz,ix) = num;
            num = num+1;
        end
    end
end

for ix = 1:size(Number_Phase,2)
    for iz = 1:size(Number_Phase,1)
        if Number_Phase(iz,ix) == 3
            Number_ind(iz,ix) = num;
            num = num+1;
        end
    end
end

num_eqns = num-1;
ind_Vx =     find(Number_Phase==1);   Number_Vx(find(Number_Vx==0))  =  Number_ind(ind_Vx);
ind_Vz =     find(Number_Phase==2);   Number_Vz(find(Number_Vz==0))  =  Number_ind(ind_Vz);
ind_P  =     find(Number_Phase==3);   Number_P (find(Number_P ==0))  =  Number_ind(ind_P );



for ix =1:2:xnum+xnum-1, for iz = 1 :2:znum+znum-1 , Number_Phase(iz,ix) = 0; end;end
tnum = 1;
for ix = 1:size(Number_Phase,2)
    for iz = 1:size(Number_Phase,1)
        if Number_Phase(iz,ix) == 0
            Number_ind(iz,ix)  = tnum;
            tnum = tnum+1;
        end
    end
end

tnum_eqns = tnum-1;
ind_T =     find(Number_Phase==0); Number_T(find(Number_T == 0))   = Number_ind(ind_T);
%  Setup  the  stiffness  Add_coeffs
A =     sparse(num_eqns,num_eqns);
Rhs_vec =     zeros(num_eqns,1);
gamma = -1e40;% punish factor
%  Setup  the  incompressibility  equations------------------------------------

%% Conservation 
ind_list  = [];
ind_val = [];
PStokes_eqn = Number_P(:);
%dVx/dx
[ind_list,ind_val] = Add_coeffs(ind_list,ind_val,Number_Vx(:,2:end   ),  ( 1/dx));
[ind_list,ind_val] = Add_coeffs(ind_list,ind_val,Number_Vx(:,1:end-1 ),  (-1/dx));
    
%dVz/dz
[ind_list,ind_val] = Add_coeffs(ind_list,ind_val,Number_Vz(2:end,:   ),  ( 1/dz));
[ind_list,ind_val] = Add_coeffs(ind_list,ind_val,Number_Vz(1:end-1,: ),  (-1/dz));

% -P/gamma
[ind_list,ind_val] = Add_coeffs(ind_list,ind_val,Number_P (:,:),  (-1/gamma));
for  i=1:size(ind_list,2)
     A  =  A  +  sparse(PStokes_eqn,ind_list(:,i),ind_val(:,i),num_eqns,num_eqns);
end



%% Vx-stokes
ind_list = [];
ind_val= [];
VxStokes_eqn = Number_Vx(2:end-1,2:end-1);
%d2Vx/dx2
[ind_list,ind_val] = Add_coeffs(ind_list , ind_val, Number_Vx(2:end-1,3:end  ),    ( 2* mu_n(2:end-1,2:end  )/dx/dx));
[ind_list,ind_val] = Add_coeffs(ind_list , ind_val, Number_Vx(2:end-1,2:end-1),    (-2*(mu_n(2:end-1,2:end)+mu_n(2:end-1,1:end-1))/dx/dx));
[ind_list,ind_val] = Add_coeffs(ind_list , ind_val, Number_Vx(2:end-1,1:end-2),    ( 2* mu_n(2:end-1,1:end-1)/dx/dx));
%d2Vx/dz2
[ind_list,ind_val] = Add_coeffs(ind_list , ind_val, Number_Vx(3:end  ,2:end-1),    (   mu_s(3:end-1,2:end-1)/dz/dz));
[ind_list,ind_val] = Add_coeffs(ind_list , ind_val, Number_Vx(2:end-1,2:end-1),    (- (mu_s(2:end-2,2:end-1)+mu_s(3:end-1,2:end-1))/dz/dz));
[ind_list,ind_val] = Add_coeffs(ind_list , ind_val, Number_Vx(1:end-2,2:end-1),    (   mu_s(2:end-2,2:end-1)/dz/dz));
%d2Vz/dx/dz
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vz(3:end-1,2:end  ),    ( mu_s(3:end-1,2:end-1)/dx/dz));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vz(2:end-2,2:end  ),    (-mu_s(2:end-2,2:end-1)/dx/dz));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vz(3:end-1,1:end-1),    (-mu_s(3:end-1,2:end-1)/dx/dz));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vz(2:end-2,1:end-1),    ( mu_s(2:end-2,2:end-1)/dx/dz));
%-dP/dx
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_P(2:end-1,2:end  ), (-1/dx));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_P(2:end-1,1:end-1), ( 1/dx));
%  Add  local  equations  to  global  Add_coeffs
for  i=1:size(ind_list,2)
    A  =  A  +  sparse(VxStokes_eqn,ind_list(:,i),ind_val(:,i),num_eqns,num_eqns);
end
% figure(1);spy(A);

%% Vz-stokes
ind_list = [];
ind_val= [];
VzStokes_eqn = Number_Vz(2:end-1,2:end-1);
%d2Vz/dz2
[ind_list,ind_val] = Add_coeffs(ind_list , ind_val, Number_Vz(3:end  ,2:end-1),    ( 2* mu_n(2:end  ,2:end-1)/dz/dz));
[ind_list,ind_val] = Add_coeffs(ind_list , ind_val, Number_Vz(2:end-1,2:end-1),    (-2*(mu_n(2:end  ,2:end-1)+mu_n(1:end-1,2:end-1))/dz/dz));
[ind_list,ind_val] = Add_coeffs(ind_list , ind_val, Number_Vz(1:end-2,2:end-1),    ( 2* mu_n(1:end-1,2:end-1)/dz/dz));
%d2Vz/dx2
[ind_list,ind_val] = Add_coeffs(ind_list , ind_val, Number_Vz(2:end-1,3:end  ),    (    mu_s(2:end-1,3:end-1)/dx/dx));
[ind_list,ind_val] = Add_coeffs(ind_list , ind_val, Number_Vz(2:end-1,2:end-1),    (-  (mu_s(2:end-1,3:end-1)+mu_s(2:end-1,2:end-2))/dx/dx));
[ind_list,ind_val] = Add_coeffs(ind_list , ind_val, Number_Vz(2:end-1,1:end-2),    (    mu_s(2:end-1,2:end-2)/dx/dx));
%d2Vz/dx/dz
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vx(2:end  ,3:end-1),    ( mu_s(2:end-1,3:end-1)/dx/dz));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vx(1:end-1,3:end-1),    (-mu_s(2:end-1,3:end-1)/dx/dz));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vx(2:end  ,2:end-2),    (-mu_s(2:end-1,2:end-2)/dx/dz));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vx(1:end-1,2:end-2),    ( mu_s(2:end-1,2:end-2)/dx/dz));
%-dP/dx
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_P(2:end  ,2:end-1), (-1/dz));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_P(1:end-1,2:end-1), ( 1/dz));
%  Add  local  equations  to  global  Add_coeffs
for  i=1:size(ind_list,2)
    A  =  A  +  sparse(VzStokes_eqn,ind_list(:,i),ind_val(:,i),num_eqns,num_eqns);
end

%% Boundary Condition
% Vx-Stokes:top
ind_list  =  [];
ind_val   =  [];
VxStokes_eqn=Number_Vx(1,2:end-1);
%d2Vx/dx2
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vx(1,3:end  ),     ( 2* mu_n(1,2:end  )/dx/dx));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vx(1,2:end-1),     (-2*(mu_n(1,1:end-1)+mu_n(1,2:end))/dx/dx));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vx(1,1:end-2),     ( 2* mu_n(1,1:end-1)/dx/dx));
%d2Vx/dz2
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vx(2,2:end-1),     (   mu_s(2,2:end-1)/dz/dz));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vx(1,2:end-1),     (  -mu_s(2,2:end-1)/dz/dz));
%d2Vz/dx/dz
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vz(2,2:end  ),     (   mu_s(2,2:end-1)/dx/dz));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vz(1,2:end  ),     (  -mu_s(2,2:end-1)/dx/dz));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vz(2,1:end-1),     (  -mu_s(1,2:end-1)/dx/dz));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vz(1,1:end-1),     (   mu_s(1,2:end-1)/dx/dz));
%-dP/dx
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_P(1,2:end  ),     (-1/dx));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_P(1,1:end-1),     ( 1/dx));
%  Add  local  equations  to  global  Add_coeffs
for  i=1:size(ind_list,2)
    A  =  A  +  sparse(VxStokes_eqn,ind_list(:,i),ind_val(:,i),num_eqns,num_eqns);
end
% Vx-Stokes:bottom
ind_list  =  [];
ind_val   =  [];
VxStokes_eqn=Number_Vx(end,2:end-1);
%d2Vx/dx2
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vx(end,3:end   ),     ( 2* mu_n(end,2:end   )/dx/dx));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vx(end,2:end-1 ),     (-2*(mu_n(end,1:end-1 )+mu_n(end,2:end))/dx/dx));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vx(end,1:end-2 ),     ( 2* mu_n(end,1:end-1 )/dx/dx));
%d2Vx/dz2
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vx(end,  2:end-1),    (  -mu_s(end-1,2:end-1)/dz/dz));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vx(end-1,2:end-1),    (   mu_s(end-1,2:end-1)/dz/dz));
%d2Vz/dx/dz
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vz(end,  2:end  ),    (   mu_s(end,  2:end-1)/dx/dz));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vz(end-1,2:end  ),    (  -mu_s(end,  2:end-1)/dx/dz));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vz(end,  1:end-1),    (  -mu_s(end-1,2:end-1)/dx/dz));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vz(end-1,1:end-1),    (   mu_s(end-1,2:end-1)/dx/dz));
%-dP/dx
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_P(end,2:end  ),     (-1/dx));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_P(end,1:end-1),     ( 1/dx));
%  Add  local  equations  to  global  Add_coeffs
for  i=1:size(ind_list,2)
    A  =  A  +  sparse(VxStokes_eqn,ind_list(:,i),ind_val(:,i),num_eqns,num_eqns);
end


% Vz-Stokes:left
ind_list  =  [];
ind_val   =  [];
VzStokes_eqn=Number_Vz(2:end-1,1);
%d2Vz/dz2
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vz(3:end,  1 ),     ( 2* mu_n(2:end,  1 )/dz/dz));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vz(2:end-1,1 ),     (-2*(mu_n(1:end-1,1 )+mu_n(2:end,1))/dz/dz));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vz(1:end-2,1 ),     ( 2* mu_n(1:end-1,1 )/dz/dz));
%d2Vz/dx2
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vz(2:end-1, 2),     (   mu_s(2:end-1, 2)/dx/dx));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vz(2:end-1, 1),     (  -mu_s(2:end-1, 2)/dx/dx));
%d2Vx/dx/dz
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vx(2:end,  2),      (   mu_s(2:end-1,  2)/dx/dz));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vx(1:end-1,2),      (  -mu_s(2:end-1,  2)/dx/dz));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vx(2:end,  1),      (  -mu_s(2:end-1,  1)/dx/dz));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vx(1:end-1,1),      (   mu_s(2:end-1,  1)/dx/dz));
% dP/dz
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_P(2:end,1  ),   (-1/dz));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_P(1:end-1,1),   ( 1/dz));
%  Add  local  equations  to  global  Add_coeffs
for  i=1:size(ind_list,2)
    A  =  A  +  sparse(VzStokes_eqn,ind_list(:,i),ind_val(:,i),num_eqns,num_eqns);
end
% Vz-Stokes:right
ind_list  =  [];
ind_val   =  [];
VzStokes_eqn=Number_Vz(2:end-1,end);
%d2Vz/dz2
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vz(3:end,  end  ),    ( 2* mu_n(2:end,  end  )/dz/dz));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vz(2:end-1,end  ),    (-2*(mu_n(1:end-1,end  )+mu_n(2:end,1))/dz/dz));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vz(1:end-2,end  ),    ( 2* mu_n(1:end-1,end  )/dz/dz));
%d2Vz/dx2
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vz(2:end-1,end  ),    (  -mu_s(2:end-1,end  )/dx/dx));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vz(2:end-1,end-1),    (   mu_s(2:end-1,end  )/dx/dx));
%d2Vx/dx/dz
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vx(2:end,  end  ),     (  mu_s(3:end,  end  )/dx/dz));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vx(1:end-1,end  ),     ( -mu_s(3:end,  end  )/dx/dz));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vx(2:end,  end-1),     ( -mu_s(2:end-1,end-1)/dx/dz));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_Vx(1:end-1,end-1),     (  mu_s(2:end-1,end-1)/dx/dz));
% dP/dz
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_P(2:end,end  ),   (-1/dz));
[ind_list,ind_val]  =  Add_coeffs(ind_list,ind_val,  Number_P(1:end-1,end),   ( 1/dz));
%  Add  local  equations  to  global  Add_coeffs
for  i=1:size(ind_list,2)
    A  =  A  +  sparse(VzStokes_eqn,ind_list(:,i),ind_val(:,i),num_eqns,num_eqns);
end



%% Rhight hand side vector
Gravity         =       9.8;
Rhs_vec(Number_Vz) =    Rhos.*Gravity;


%% Apply boundary conditions
% bottom boundary:Vz=0
bcind_bot           =       Number_Vz(end,:);
bcval_bot           =       0*ones(size(bcind_bot));
% top boundary:Vz=0
bcind_top           =       Number_Vz(1,:);
bcval_top           =       0*ones(size(bcind_top));
% left boundary:Vx=0
bcind_left          =       Number_Vx(:,1);
bcval_left          =       0*ones(size(bcind_left));
% right boundary:Vx=0
bcind_right         =       Number_Vx(:,end);
bcval_right         =       0*ones(size(bcind_right));

Bc_ind          =       [bcind_left',bcind_bot,bcind_right',bcind_top];
Bc_val          =       [bcval_left',bcval_bot,bcval_right',bcval_top];

% A(sub2ind(size(A),Bc_ind(:),Bc_ind(:))) =   1;
% Rhs_vec(Bc_ind,1)   =       Bc_val;
% the below is to make the A matrice symmetric: move the known values to
% right side
A(Bc_ind,:)         =       0;
Rhs_vec             =       Rhs_vec-A(:,Bc_ind)*Bc_val';
A(:,Bc_ind)         =       0;
store_f             =       spdiags(A,0);
store_f(Bc_ind)     =       1;
A                   =       spdiags(store_f,0,A);
Rhs_vec(Bc_ind)     =       Bc_val;


%% Testing if A is symmetric
figure(11),clf
spy(A)
diffA=A-A';
max(max(abs(diffA)))
% print(gcf,'-djpeg','-r300','-opengl','valvis52.jpg')

%% Solving equations
C               =       A\Rhs_vec;
Vx              =       C(Number_Vx);
Vz              =       C(Number_Vz);
Pressure        =       C(Number_P);

SecYear         =       3600*24*365.25;
cmYear          =       SecYear*1e2;
Vx              =       reshape(Vx,size(XVx,1),size(XVx,2)).*cmYear;
Vz              =       reshape(Vz,size(XVz,1),size(XVz,2)).*cmYear;
Pressure        =       reshape(Pressure,size(XP,1),size(XP,2));

%% recalculate viscosity
etaxx         =       (Vx(:,1:end-1)-Vx(:,2:end))./dx./cmYear;
etayy         =       (Vz(1:end-1,:)-Vz(2:end,:))./dz/cmYear;
etaxy         =       1/2.*((Vx(:,1:end-1)-Vx(:,2:end))./dz+(Vz(1:end-1,:)-Vz(2:end,:))./dx)/cmYear;
etaII         =       sqrt(1/2*(etaxx.^2+etayy.^2)+etaxy.^2);
e1            =       ((Ea+Va.*Pressure)./(R.*Tn));
% mu_n          =       F2.*mu_n.*1./(h^m.*etaII.^(n-1)).*exp(e1);

%% Temperature equation
B  = sparse(tnum_eqns,tnum_eqns);
TRhs_vec = zeros(tnum_eqns,1);
dtmax = (Cp(1)*rho(1))./k(1)/(2*(1/dx^2+1/dz^2));
dt = (Cp(1)*rho(1))./k(1)/(2*(1/dx^2+1/dz^2));
% dt = 1;
k_x   = 2*(k_n(:,1:end-1).*k_n(:,2:end))./(k_n(:,1:end-1)+k_n(:,2:end));
k_z   = 2*(k_n(1:end-1,:).*k_n(2:end,:))./(k_n(1:end-1,:)+k_n(2:end,:));

ind_list = [];
ind_val  = [];
Temp_eqns = Number_T(2:end-1,2:end-1);
% Tn5
[ind_list,ind_val]  = Add_coeffs(ind_list, ind_val, Number_T(2:end-1,3:end), (-k_z(2:end,2:end-1)/dx/dx));
% Tn3
[ind_list,ind_val]  = Add_coeffs(ind_list, ind_val, Number_T(2:end-1,2:end-1), (Rho(2:end-1,2:end-1).*Cp_n(2:end-1,2:end-1)/dt));
[ind_list,ind_val]  = Add_coeffs(ind_list, ind_val, Number_T(2:end-1,2:end-1), ((k_x(2:end-1,1:end-1)+k_x(2:end-1,2:end))/dz/dz+(k_z(1:end-1,2:end-1)+k_z(2:end,2:end-1))/dx/dx));
% Tn1
[ind_list,ind_val]  = Add_coeffs(ind_list, ind_val, Number_T(2:end-1,1:end-2), (-k_z(2:end  ,2:end-1)/dx/dx));
% Tn4
[ind_list,ind_val]  = Add_coeffs(ind_list, ind_val, Number_T(3:end  ,2:end-1), (-k_x(2:end-1,2:end  )/dz/dz));
% Tn2
[ind_list,ind_val]  = Add_coeffs(ind_list, ind_val, Number_T(1:end-2,2:end-1), (-k_x(2:end-1,1:end-1)/dz/dz));
% T03
% [ind_list,ind_val]  = Add_coeffs(ind_list, ind_val, Number_T(2:end-1,2:end-2), (Rho(2:end-1,2:end-1).*Cp_n(2:end-1,2:end-2)/dtmax));

for  i=1:size(ind_list,2)
    B  =  B  +  sparse(Temp_eqns,ind_list(:,i),ind_val(:,i),tnum_eqns,tnum_eqns);
end
% figure(12),clf
% spy(B)

%% Temperature Boundary condition
ind_list = [];
ind_val  = [];
Temp_eqns= Number_T(:,1);
% Temperature left : An insulating boundary condition
[ind_list,ind_val]  = Add_coeffs(ind_list, ind_val, Number_T(:,1), (-k_x(:,1)/dx));
[ind_list,ind_val]  = Add_coeffs(ind_list, ind_val, Number_T(:,2), ( k_x(:,1)/dx));

for  i=1:size(ind_list,2)
    B  =  B  +  sparse(Temp_eqns,ind_list(:,i),ind_val(:,i),tnum_eqns,tnum_eqns);
end

ind_list = [];
ind_val  = [];
Temp_eqns= Number_T(:,end);
% Temperature right : An insulating boundary condition
[ind_list,ind_val]  = Add_coeffs(ind_list, ind_val, Number_T(:,end)  , ( k_x(:,end)/dx));
[ind_list,ind_val]  = Add_coeffs(ind_list, ind_val, Number_T(:,end-1), (-k_x(:,end)/dx));
for  i=1:size(ind_list,2)
    B  =  B  +  sparse(Temp_eqns,ind_list(:,i),ind_val(:,i),tnum_eqns,tnum_eqns);
end
Tnold  = Ts;
% TRhs_vec(Number_T) =    Tnold(2:end-1,2:end-1).*Rho(2:end-1,2:end-1).*Cp_n(2:end-1,2:end-1)/dt;
TRhs_vec(Number_T) =    Tnold.*Rho.*Cp_n/dt;

%% Apply Temperature Boundary condition
% Temperature top : constant temperature
bcind_top           = Number_T(1,:);
bcval_top           = Ttop*ones(size(bcind_top));
% Temperature bottom: constant temperature
bcind_bot           = Number_T(end,:);
bcval_bot           = Tbot*ones(size(bcind_bot));

Bc_ind          =       [bcind_bot,bcind_top];
Bc_val          =       [bcval_bot,bcval_top];

B(Bc_ind,:)         =       0;
TRhs_vec            =       TRhs_vec-B(:,Bc_ind)*Bc_val';
B(:,Bc_ind)         =       0;
store_f             =       spdiags(B,0);
store_f(Bc_ind)     =       1;
B                   =       spdiags(store_f,0,B);
TRhs_vec(Bc_ind)    =       Bc_val;

D                   =       B\TRhs_vec;
T                   =       TRhs_vec(Number_T);

% Tnew                =       T;
Ts                  =       T;

%% Plotting
figure(1),clf
subplot(121)
pcolor(XVx,ZVx,Vx)
shading interp;
colorbar
title('Vx');
axis equal tight
subplot(122)
pcolor(XVz,ZVz,Vz)
shading interp;
colorbar
title('Vz')  
axis equal tight
% print(gcf,'-djpeg','-r300','-opengl','valvis1.jpg')

Vx_cen = (Vx(:,1:end-1)+Vx(:,2:end))/2;
Vz_cen = (Vz(1:end-1,:)+Vz(2:end,:))/2;

figure(2),clf
pcolor(XP,ZP,Pressure)
shading interp
colorbar
hold on
quiver(XP,ZP,Vx_cen,Vz_cen,'k')
title('Pressure, Velocity(arrow)')  
axis equal tight
% print(gcf,'-djpeg','-r300','-opengl','valvis2.jpg')

figure(3),clf
rho=(Rhos(1:end-1,:)+Rhos(2:end,:))/2.0;
pcolor(XP,ZP,rho);
shading interp
colorbar
hold on
quiver(XP,ZP,Vx_cen,Vz_cen,'k')
title('Density, Velocity(arrow)') 
axis equal tight
% print(gcf,'-djpeg','-r300','-opengl','valvis3.jpg')

figure(4),clf
nu = mu_n(:,:);
pcolor(XP,ZP,nu);
shading interp
colorbar
hold on
quiver(XP,ZP,Vx_cen,Vz_cen,'k')
title('Viscosity, Velocity(arrow)') 
axis equal tight

% print(gcf,'-djpeg','-r300','-opengl','valvis4.jpg')

figure(5),clf
pcolor(xgrid,zgrid,Ts)
shading interp
colorbar
hold on
% quiver(Px,Pz,Vx_cen,Vz_cen,'k')
title('Temperature')  
axis equal tight
% print(gcf,'-djpeg','-r300','-opengl','valvis6.jpg')

figure(6),clf
pcolor(XP,ZP,etaII)
shading interp
colorbar
hold on
% quiver(Px,Pz,Vx_cen,Vz_cen,'k')
title('\epsilonII')  
axis equal tight
% print(gcf,'-djpeg','-r300','-opengl','valvis7.jpg')


toc;