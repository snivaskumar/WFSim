function [StrucDiscretization,StrucBCs] = BoundaryConditions(Nx,Ny,StrucDiscretization,u,v,Linearversion)
ax = StrucDiscretization.ax;
ay = StrucDiscretization.ay;

bbx = sparse((Nx-3)*(Ny-2),(Nx-3)*(Ny-2)+(Nx-2)*(Ny-3)+(Nx-2)*(Ny-2));
bby = sparse((Nx-2)*(Ny-3),(Nx-3)*(Ny-2)+(Nx-2)*(Ny-3)+(Nx-2)*(Ny-2));

% Zero gradient outflow
% for u-direction
ax.aP((Nx-1),(2:Ny-1))  = ax.aP((Nx-1),(2:Ny-1)) - ax.aE((Nx-1),(2:Ny-1)); %NORTH
ax.aP((1:Nx-1),(Ny-1))  = ax.aP((1:Nx-1),(Ny-1)) - ax.aN((1:Nx-1),(Ny-1)); %EAST
ax.aP((1:Nx-1),(2))     = ax.aP((1:Nx-1),2)      - ax.aS((1:Nx-1),2);

%% for v-direction
ay.aP((Nx-1),(1:Ny))  = ay.aP((Nx-1),(1:Ny)) - ay.aE((Nx-1),(1:Ny));
ay.aP((1:Nx),(Ny-1))  = ay.aP((1:Nx),(Ny-1)) - ay.aN((1:Nx),(Ny-1));
ay.aP((1:Nx),(3))     = ay.aP((1:Nx),3)      - ay.aS((1:Nx),3); % changed to 3 3 2 instead of 2 2 1

% Inflow boundary for non linear model
bx      = kron([1;zeros(Nx-4,1)],(ax.aW(3,2:end-1).*u(2,2:end-1))'); %changed to 3: 2 instead of 2:2
by      = [v(1,3:Ny-1)'.*ay.aW(2,3:Ny-1)';zeros((Nx-3)*(Ny-3),1)]; %changed to 2:3 inst

% Write nonlinear to outputs
StrucDiscretization.ax  = ax;
StrucDiscretization.ay  = ay;
StrucBCs.bx             = bx;
StrucBCs.by             = by;
end