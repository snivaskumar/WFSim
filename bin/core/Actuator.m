function [output,Ueffect,a,Power,CT] = Actuator(Wp,input,sol,options)

Nx          = Wp.mesh.Nx;
Ny          = Wp.mesh.Ny;
dyy2        = Wp.mesh.dyy2;
xline       = Wp.mesh.xline;
yline       = Wp.mesh.yline;
ylinev      = Wp.mesh.ylinev;
Rho         = Wp.site.Rho;
Drotor      = Wp.turbine.Drotor;
powerscale  = Wp.turbine.powerscale;
N           = Wp.turbine.N;
F           = Wp.turbine.forcescale; % http://www.nrel.gov/docs/fy05osti/36834.pdf
Ar          = pi*(0.5*Drotor)^2;
scale       = 3;                                  % To scale the force in the y-direction
k           = options.k;

[Sm.x,Sm.dx]    = deal(sparse(Nx-3,Ny-2));            % Input x-mom nonlinear and linear
[Sm.y,Sm.dy]    = deal(sparse(Nx-2,Ny-3));            % Input y-mom nonlinear and linear
[Sm.xx,Sm.dxx]  = deal(sparse((Nx-3)*(Ny-2),2*N));    % Input x-mom nonlinear and linear qlpv
[Sm.yy,Sm.dyy]  = deal(sparse((Nx-2)*(Ny-3),2*N));    % Input y-mom nonlinear and linear qlpv


for kk=1:N
    % Clear for multiple turbine case
    tempx = sparse(Nx-3,Ny-2);
    tempy = sparse(Nx-2,Ny-3);
    
    x  = xline(kk,:);  % Turbine x-pos in field
    y  = yline{kk};    % Turbine y-pos in field
    yv = ylinev{kk};   % Corrected turbine y-pos in field
    
    %%
    vv            = 0.5*diff(sol.v(x,yv))+sol.v(x,yv(1:end-1)); % Bart: this can be fixed!
    uu            = sol.u(x,y);
    U{kk}         = sqrt(uu.^2+vv.^2);
    phi{kk}       = atan(vv./uu);
    Ue{kk}        = cos(phi{kk}+input.phi(kk)/180*pi).*U{kk};
  
    a(kk)         = input.beta(kk)/(input.beta(kk)+1);  
    CT(kk)        = 4*a(kk)*F*(1-a(kk));
   
    Ueffect(kk)   = mean(Ue{kk})/(1-a(kk));     % Estimation effective wind speed
    
    %% Thrust force
    Fthrust       = 1/2*Rho*Ue{kk}.^2*CT(kk)*(input.beta(kk)+1).^2;
    Fx            = Fthrust.*cos(input.phi(kk)*pi/180);
    Fy            = Fthrust.*sin(input.phi(kk)*pi/180);
    
    %% Power
    Power(kk)     = mean(powerscale*.5*Rho*Ar*(Ue{kk}).^3*CT(kk)*cos(input.phi(kk)*pi/180)^(1.88));

    if k==2
        a0      = a(kk);            % linearisation point
        Uinf0   = Ueffect(kk);      % linearisation point
        
        CreateTurbineModel
    end    
    
    if k>=2
        w         = [ zeros(1,1) ; input.dbeta(kk)];
        
        xT1(:,k+1) = A*xT1(:,k) + B*w(:,k);
        
        x(:,k+1) = A*x(:,k) + B*w(:,k);
        y(:,k)   = C*x(:,k) + D*w(:,k);
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    %% Input to Ax=b
    Sm.x(x-2,y-1)           = -Fx'.*dyy2(1,y)';                                                                  % Input x-mom nonlinear                           % Input x-mom linear
    Sm.y(x-1,y(2:end)-2)    = scale*Fy(2:end)'.*dyy2(1,y(2:end))';                                               % Input y-mom nonlinear
    
end
% Write to output
output.Sm  = Sm;


