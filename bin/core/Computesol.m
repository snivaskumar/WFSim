function [sol,sys] = Computesol(sys,input,sol,k,it,options)

sol.x(sys.pRCM,1) = sys.A(sys.pRCM,sys.pRCM)\sys.b(sys.pRCM);

end