function [output] = SpatialDiscr_Hybrid(Wp,u,v,Linearversion)
% dxx   = \Delta x_{I,I+1}
% dyy   = \Delta y_{J,J+1}
% dxx2  = \Delta x_{i,i+1}
% dyy2  = \Delta y_{j,j+1}
% ldxx  = I
% ldyy  = J
% ldxx2 = i
% ldyy2 = j

Nx     = Wp.mesh.Nx;
Ny     = Wp.mesh.Ny;
dxx    = Wp.mesh.dxx;
dyy    = Wp.mesh.dyy;
dxx2   = Wp.mesh.dxx2;
dyy2   = Wp.mesh.dyy2;
xline  = Wp.mesh.xline;
yline  = Wp.mesh.yline;

Rho    = Wp.site.Rho;
mu     = Wp.site.mu;
lmu    = Wp.site.lmu;
Turb   = Wp.site.turbul;
m      = Wp.site.m;
n      = Wp.site.n;

Drotor = Wp.turbine.Drotor;
N      = Wp.turbine.N;

% Init
[ax.aE,ax.aW,ax.aS,ax.aN,ax.aP]         = deal(zeros(Nx,Ny));
[Fex,Fwx,Fsx,Fnx,dFex,dFwx,dFnx,dFsx]   = deal(zeros(Nx,Ny));
[ay.aE,ax.aW,ay.aS,ay.aN,ay.aP]         = deal(zeros(Nx,Ny));
[Fey,Fwy,Fsy,Fny,dFey,dFwy,dFny,dFsy]   = deal(zeros(Nx,Ny));

%%  Setting the coefficients according to the hybrid differencing scheme %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% x-direction
% Define diffusion coefficients
Dxe             = mu./[dxx2(2:end,:);dxx2(end,:)].*dyy2;
Dxw             = mu./dxx2.*dyy2;
Dxn             = mu./[dyy(:,2:end) dyy(:,end)].*dxx;
Dxs             = mu./dyy.*dxx;

% Define convection coefficients and its derivatives
% Fex = c ( u_{i,J} + u_{i+1,J} )
dFex(1:Nx-1,1:Ny) = Rho*0.5*dyy2(1:Nx-1,1:Ny);
Fex(1:Nx-1,1:Ny)  = dFex(1:Nx-1,1:Ny).*( u(2:Nx,1:Ny) + u(1:Nx-1,1:Ny) );
% Few = c ( u_{i,J} + u_{i-1,J} )
dFwx(2:Nx,1:Ny)   = Rho*0.5*dyy2(2:Nx,1:Ny);
Fwx(2:Nx,1:Ny)    = dFwx(2:Nx,1:Ny).*( u(2:Nx,1:Ny) + u(1:Nx-1,1:Ny) ); %Zelfde als Fex?
% Fnx = c ( v_{I-1,j+1} + v_{I,j+1} )
dFnx(2:Nx,1:Ny-1) = Rho*0.5*dxx(2:Nx,1:Ny-1);
Fnx(2:Nx,1:Ny-1)  = dFnx(2:Nx,1:Ny-1).*( v(2:Nx,2:Ny) + v(1:Nx-1,2:Ny) );
% Fsx = c ( v_{I-1,j} + v_{I,j} )
dFsx(2:Nx,1:Ny)   = Rho*0.5*dxx(2:Nx,1:Ny);
Fsx(2:Nx,1:Ny)    = dFsx(2:Nx,1:Ny).*( v(2:Nx,1:Ny) + v(1:Nx-1,1:Ny)); % Waarom deze een andere size dan de andere drie?

ax.aE             = max(max(-Fex,Dxe-0.5*Fex),zeros(Nx,Ny));
ax.aW             = max(max(Fwx,Dxw+0.5.*Fwx),zeros(Nx,Ny));
ax.aN             = max(max(-Fnx,Dxn-0.5*Fnx),zeros(Nx,Ny));
ax.aS             = max(max(Fsx,Dxs+0.5*Fsx),zeros(Nx,Ny));
ax.aP             = ax.aW + ax.aE + ax.aS + ax.aN + Fex - Fwx + Fnx - Fsx;

%% y-direction
% Define diffusion coefficients
Dye             = mu./[dxx(2:end,:);dxx(end,:)].*dyy;
Dyw             = mu./dxx.*dyy;
Dyn             = mu./[dyy2(:,2:end) dyy2(:,end)].*dxx2;
Dys             = mu./dyy2.*dxx2;

% Define convection coefficients and its derivatives
% Fey = c ( u_{i+1,J} + u_{i+1,J-1} )
dFey(1:Nx-1,2:Ny) = Rho*0.5*dyy(1:Nx-1,2:Ny);
Fey(1:Nx-1,2:Ny)  = dFey(1:Nx-1,2:Ny).*( u(2:Nx,2:Ny) + u(2:Nx,1:Ny-1) );
% Fwy = c ( u_{i,J} + u_{i,J-1} )
dFwy(1:Nx,2:Ny)   = Rho*0.5*dyy(1:Nx,2:Ny);
Fwy(1:Nx,2:Ny)    = dFwy(1:Nx,2:Ny).*( u(1:Nx,2:Ny) + u(1:Nx,1:Ny-1) );
% Fny = c ( v_{I,j+1} + v_{I,j} )
dFny(1:Nx,1:Ny-1) = Rho*0.5*dxx2(1:Nx,1:Ny-1);
Fny(1:Nx,1:Ny-1)  = dFny(1:Nx,1:Ny-1).*( v(1:Nx,1:Ny-1) + v(1:Nx,2:Ny) );
% Fsy = c ( v_{I,j-1} + v_{I,j} )
dFsy(1:Nx,2:Ny)   = Rho*0.5*dxx2(1:Nx,2:Ny);
Fsy(1:Nx,2:Ny)    = dFsy(1:Nx,2:Ny).*( v(1:Nx,1:Ny-1) + v(1:Nx,2:Ny) );

ay.aE             = max(max(-Fey,Dye-0.5*Fey),zeros(Nx,Ny));
ay.aW             = max(max(Fwy,Dyw+0.5.*Fwy),zeros(Nx,Ny));
ay.aN             = max(max(-Fny,Dyn-0.5*Fny),zeros(Nx,Ny));
ay.aS             = max(max(Fsy,Dys+0.5*Fsy),zeros(Nx,Ny));
ay.aP             = ay.aW + ay.aE + ay.aS + ay.aN + Fey - Fwy + Fny - Fsy;

% Turbulence model
if Turb==1    
    Turbulence    
end

output.ax = ax;
output.ay = ay;