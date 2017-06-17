% Linear turbine model
clear;clc;close all;
Rho  = 1.2;
Ar   = pi*45^2;
cf   = 2;
Uinf = 8;%linspace(0,10,50);
a    = linspace(0,.7,50);

for kk=1:numel(Uinf)   
    for ll=1:numel(a)
 
        CT(kk,ll)       =  4*a(ll)*cf*(1-a(ll));
        CP(kk,ll)       =  4*a(ll)*cf*(1-a(ll))^2;        
        F(kk,ll)        = CT(kk,ll)*.5*Rho*Ar*( Uinf(kk) )^2;        
        P(kk,ll)        = CP(kk,ll)*.5*Rho*Ar*( Uinf(kk) )^3;
       
    end
end

figure(1);clf;
subplot(1,2,1)
plot(a,F);grid
xlabel('a');ylabel('F');
subplot(1,2,2)
plot(a,P);grid
xlabel('a');ylabel('P');
%% Steady-state linear model
clear;clc;close all;
Rho  = 1.2;
Ar   = pi*45^2;
cf   = 2;
Uinf = linspace(3,8,50);
a    = .3*ones(1,50);

a0      = .3;
Uinf0   = 6;
CT0     = 4*a0*cf*(1-a0);
CP0     = 4*a0*cf*(1-a0)^2;
F0      = CT0*.5*Rho*Ar*( Uinf0 )^2;
P0      = CP0*.5*Rho*Ar*( Uinf0 )^3;
dFdU    = CT0*.5*Rho*Ar*2*( Uinf0 );
dCTda   = 4*cf-8*cf*a0;
dFda    = dCTda*.5*Rho*Ar*( Uinf0 )^2;
dPdU    = CP0*.5*Rho*Ar*3*( Uinf0 )^2;
dCPda   = 4*cf*(1-a0)^2-8*a0*cf*(1-a0);
dPda    = dCPda*.5*Rho*Ar*( Uinf0 )^3;

da      = 0*linspace(-.2,.2,50);
dU      = linspace(-1,1,50);

for kk=1:numel(da)
    
    CT(kk) = 4*a(kk)*cf*(1-a(kk));
    CP(kk) = 4*a(kk)*cf*(1-a(kk))^2;   
    F(kk)  = CT(kk)*.5*Rho*Ar*( Uinf(kk) )^2;
    P(kk)  = CP(kk)*.5*Rho*Ar*( Uinf(kk) )^3;   
    dF(kk) = dFdU*dU(kk) + dFda*da(kk);
    dP(kk) = dPdU*dU(kk) + dPda*da(kk);
end

figure(2);clf;
subplot(1,2,1)
plot(Uinf,F);hold on;
plot(Uinf0+dU,F0+dF);
xlabel('U\infty');ylabel('F');grid;
subplot(1,2,2)
plot(Uinf,P);hold on;
plot(Uinf0+dU,P0+dP);
xlabel('U\infty');ylabel('P');grid;

Uinf = 8*ones(1,50);
a    = linspace(0,.7,50);

a0      = .3;
Uinf0   = Uinf(1);
CT0     = 4*a0*cf*(1-a0);
CP0     = 4*a0*cf*(1-a0)^2;
F0      = CT0*.5*Rho*Ar*( Uinf0 )^2;
P0      = CP0*.5*Rho*Ar*( Uinf0 )^3;
dFdU    = CT0*.5*Rho*Ar*2*( Uinf0 );
dCTda   = 4*cf-8*cf*a0;
dFda    = dCTda*.5*Rho*Ar*( Uinf0 )^2;
dPdU    = CP0*.5*Rho*Ar*3*( Uinf0 )^2;
dCPda   = 4*cf*(1-a0)^2-8*a0*cf*(1-a0);
dPda    = dCPda*.5*Rho*Ar*( Uinf0 )^3;
da      = linspace(-.2,.2,50);
dU      = 0*linspace(-.3,.3,50);

for kk=1:numel(da)
    
    CT(kk) = 4*a(kk)*cf*(1-a(kk));
    CP(kk) = 4*a(kk)*cf*(1-a(kk))^2;   
    F(kk)  = CT(kk)*.5*Rho*Ar*( Uinf(kk) )^2;
    P(kk)  = CP(kk)*.5*Rho*Ar*( Uinf(kk) )^3;
    
    dF(kk) = dFdU*dU(kk) + dFda*da(kk);
    dP(kk) = dPdU*dU(kk) + dPda*da(kk);
end

figure(1);clf;
subplot(1,2,1)
plot(a,F);hold on;
plot(a0+da,F0+dF);
xlabel('a');ylabel('F');grid;
subplot(1,2,2)
plot(a,P);hold on;
plot(a0+da,P0+dP);
xlabel('a');ylabel('P');grid;

%% Dynamic linear model
clear;clc;close all
Rho  = 1.2;
Ar   = pi*45^2;
cf   = 2;

dU   = 0*ones(1,20);
da   = [.1*ones(1,10) -.1*ones(1,10)];

% Compute steady-state situation linear and nonlinear
a0      = .3;   % linearisation point
Uinf0   = 6;    % linearisation point
CT0     = 4*a0*cf.*(1-a0);
CP0     = 4*a0*cf.*(1-a0).^2;
F0      = CT0*.5*Rho*Ar.*( Uinf0 ).^2;
P0      = CP0*.5*Rho*Ar.*( Uinf0 ).^3;
dFdU    = CT0*.5*Rho*Ar*2.*( Uinf0 );
dCTda   = 4*cf-8*cf*a0;
dFda    = dCTda*.5*Rho*Ar.*( Uinf0 ).^2;
dPdU    = CP0*.5*Rho*Ar*3.*( Uinf0 ).^2;
dCPda   = 4*cf*(1-a0).^2-8*a0*cf.*(1-a0);
dPda    = dCPda*.5*Rho*Ar.*( Uinf0 ).^3;
dF      = dFdU.*dU + dFda.*da;
dP      = dPdU.*dU + dPda.*da;

% Transfer from w = [dU da]^T to y = [dF dP]^T
s         = tf('s');
wn        = 1;
z         = .7;
Hi        = ss(wn^2/(s^2 + 2*z*wn*s + wn^2));
H         = Hi*[dFdU dFda;dPdU dPda];
Hd        = c2d(H,1);
[A,B,C,D] = ssdata(Hd);

w  = [ dU ; da]; 
x  = zeros( size(A,1),numel(da) );

for k=1:numel(da)-1 
    x(:,k+1) = A*x(:,k) + B*w(:,k);
    y(:,k)   = C*x(:,k) + D*w(:,k);      
end
    
figure(3);clf;
subplot(2,2,1)
plot(y(1,:));hold on;
plot(dF);
xlabel('k');ylabel('dF');grid;
subplot(2,2,2)
plot(y(2,:));hold on;
plot(dP);
xlabel('k');ylabel('dP');grid;
subplot(2,2,3)
plot(w(1,:));hold on;
plot(dU);
xlabel('k');ylabel('dU');grid;
subplot(2,2,4)
plot(w(2,:));hold on;
plot(da);
xlabel('k');ylabel('da');grid;
