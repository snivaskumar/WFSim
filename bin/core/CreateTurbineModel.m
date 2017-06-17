a0      = Wp.turbine.a0{kk};
Uinf0   = Wp.turbine.Uinf0{kk};

CT0     = 4*a0*F.*(1-a0);
CP0     = 4*a0*F.*(1-a0).^2;
dFdU    = CT0*.5*Rho*2.*( Uinf0 );
dCTda   = 4*F-8*F*a0;
dFda    = dCTda*.5*Rho.*( Uinf0 ).^2;
dPdU    = CP0*.5*Rho*Ar*3.*( Uinf0 ).^2;
dCPda   = 4*F*(1-a0).^2-8*a0*F.*(1-a0);
dPda    = dCPda*.5*Rho*Ar.*( Uinf0 ).^3;

% Linarisation points output
Wp.turbine.F0{kk}   = Fx';
Wp.turbine.P0{kk}   = CP(kk)*powerscale*.5*Rho*Ar*Ueffect(kk).^3;

s         = tf('s');
wn        = .5;
z         = 2;
Hi        = ss(wn^2/(s^2 + 2*z*wn*s + wn^2));
H         = Hi*[dFdU dFda;dPdU dPda];
Hd        = c2d(H,Wp.sim.h);
[A,B,C,D] = ssdata(Hd);

% Transfer from w = [dU da]^T to y = [dF dP]^T
Wp.turbine.A{kk} = A;
Wp.turbine.B{kk} = B;
Wp.turbine.C{kk} = C;
Wp.turbine.D{kk} = D;

Wp.turbine.xT{kk} = zeros(size(A,1),Wp.sim.NN);
Wp.turbine.yT{kk} = zeros(size(C,1),Wp.sim.NN);