function [sol,eps] = MapSolution(Nx,Ny,sol,k,it,options)


sol.uu(3:end-1,2:end-1) = reshape(sol.x(1:(Nx-3)*(Ny-2)),Ny-2,Nx-3)';
sol.vv(2:end-1,3:end-1) = reshape(sol.x((Nx-3)*(Ny-2)+1:(Nx-3)*(Ny-2)+(Nx-2)*(Ny-3)),Ny-3,Nx-2)';
%pp(2:end-1,2:end-1)     = reshape([sol.x((Nx-3)*(Ny-2)+(Nx-2)*(Ny-3)+1:end);0],Ny-2,Nx-2)';
sol.pp(2:end-1,2:end-1) = reshape([sol.x((Nx-3)*(Ny-2)+(Nx-2)*(Ny-3)+1:end);0;0],Ny-2,Nx-2)';
sol.pp(isinf(sol.pp))   = 0;

% Check if solution converged
Normv{k} = norm(vec(sol.v(2:end-1,3:end-1)-sol.vv(2:end-1,3:end-1)));
Normu{k} = norm(vec(sol.u(3:end-1,2:end-1)-sol.uu(3:end-1,2:end-1)));
eps      = sqrt((Normv{k}+Normu{k}))/((Ny-2)*(Nx-2))/2;

if k==1; alpha      = min(1-.9^it,1); else alpha=1; end;
sol.u(3:end-1,2:end-1)  = (1-alpha)*sol.u(3:end-1,2:end-1)+alpha*sol.uu(3:end-1,2:end-1);
sol.v(2:end-1,3:end-1)  = (1-alpha)*sol.v(2:end-1,3:end-1)+alpha*sol.vv(2:end-1,3:end-1);
sol.p(2:end-1,2:end-1)  = (1-alpha)*sol.p(2:end-1,2:end-1)+alpha*sol.pp(2:end-1,2:end-1);

% Update velocities for next iteration and boundary conditions
[sol.u,sol.v,sol.p]    = Updateboundaries(Nx,Ny,sol.u,sol.v,sol.p);

end