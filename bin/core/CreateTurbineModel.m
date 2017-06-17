CT0     = 4*a0*F.*(1-a0);
CP0     = 4*a0*F.*(1-a0).^2;
F0      = CT0*.5*Rho*Ar.*( Uinf0 ).^2;
P0      = CP0*.5*Rho*Ar.*( Uinf0 ).^3;
dFdU    = CT0*.5*Rho*Ar*2.*( Uinf0 );
dCTda   = 4*F-8*F*a0;
dFda    = dCTda*.5*Rho*Ar.*( Uinf0 ).^2;
dPdU    = CP0*.5*Rho*Ar*3.*( Uinf0 ).^2;
dCPda   = 4*F*(1-a0).^2-8*a0*F.*(1-a0);
dPda    = dCPda*.5*Rho*Ar.*( Uinf0 ).^3;

s         = tf('s');
wn        = 1;
z         = .7;
Hi        = ss(wn^2/(s^2 + 2*z*wn*s + wn^2));
H         = Hi*[dFdU dFda;dPdU dPda];
Hd        = c2d(H,1);
[A,B,C,D] = ssdata(Hd);

Wp.turbine.A = A;
Wp.turbine.B = B;
Wp.turbine.C = C;
Wp.turbine.D = D;
