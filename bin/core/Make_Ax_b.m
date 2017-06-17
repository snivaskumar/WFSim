function [sys,Power,Ueffect,a,CT,Wp] = Make_Ax_b(Wp,sys,sol,input,B1,B2,bc,k,options)
Nx    = Wp.mesh.Nx;
Ny    = Wp.mesh.Ny;

% Decide whether to start from uniform flow field or steady state
if k == 1 && options.startUniform == 0
    dt = Inf;
else
    dt        = Wp.sim.h;
end;
options.k = k;

[StrucDiscretization]                = SpatialDiscr_Hybrid(Wp,sol.u,sol.v,0); % Spatial discretization
[StrucDiscretization,StrucDynamical] = Dynamical(Wp,StrucDiscretization,sol.u,sol.v,dt,0); % Dynamical term
[StrucActuator,Ueffect,a,Power,CT,Wp]= Actuator(Wp,input,sol,options); % Actuator
[StrucDiscretization,StrucBCs]       = BoundaryConditions(Nx,Ny,StrucDiscretization,sol.u,sol.v,0); % Zero gradient boundary conditions momentum equations

% Setup A matrix
Ay    = MakingSparseMatrix(Nx,Ny,StrucDiscretization.ay,2,3,1);
Ax    = MakingSparseMatrix(Nx,Ny,StrucDiscretization.ax,3,2,1);


sys.A = [blkdiag(Ax,Ay) [B1;B2]; [B1;2*B2]' sparse((Nx-2)*(Ny-2),(Nx-2)*(Ny-2))];

sys.M = blkdiag(...
    spdiags(StrucDynamical.ccx,0,Wp.Nu,Wp.Nu),...
    spdiags(StrucDynamical.ccy,0,Wp.Nv,Wp.Nv),...
    spdiags(zeros(Wp.Np,1),0,Wp.Np,Wp.Np));

sys.b    = [StrucBCs.bx+StrucDynamical.cx+vec(StrucActuator.Sm.x');
    StrucBCs.by+StrucDynamical.cy+vec(StrucActuator.Sm.y');
    bc];
sys.m    = [StrucBCs.bx+vec(StrucActuator.Sm.x');
    StrucBCs.by+vec(StrucActuator.Sm.y');
    bc];

sys.A(size(Ax,1)+size(Ay,1)+size(B1',1)-(Ny-2)+1,:) = [];
sys.b(size(Ax,1)+size(Ay,1)+size(B1',1)-(Ny-2)+1,:) = [];
sys.m(size(Ax,1)+size(Ay,1)+size(B1',1)-(Ny-2)+1,:) = [];
sys.A(:,size(Ax,1)+size(Ay,1)+size(B1',1)-(Ny-2)+1) = [];
sys.A(:,end) = [];sys.A(end,:) = [];sys.b(end)=[];sys.m(end)=[];

if k==1
    sys.pRCM = symrcm(sys.A);
end;
