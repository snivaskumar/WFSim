%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Types of systems when Projection=0:
% Linearversion = 0;
% sys.A * x_{k+1}  = sys.b(x_k,u_k)
% sys.A * x_{k+1}  = sys.M * x_k + sys.m(x_k,u_k)
% Linearversion = 1;
% sys.A * dx_{k+1} = sys.Al * dx_k + sys.bl(du_k)
% sys.A * dx_{k+1} = sys.Al * dx_k + sys.Bl * du_k

clear; clc; close all;

WFSim_addpaths

%% Initialize script
options.startUniform   = 0;                      % Start from a uniform flowfield (true) or a steady-state solution (false)
Wp.name                = 'RobustMpc';

Animate       = 10 ;                     % Show 2D flow fields every x iterations (0: no plots)
plotMesh      = 0;                      % Show meshing and turbine locations
conv_eps      = 1e-6;                   % Convergence threshold
max_it_dyn    = 1;                      % Maximum number of iterations for k > 1

if options.startUniform==1; max_it = 1; else max_it = 50; end

% WFSim general initialization script
[Wp,sol,sys,Power,CT,a,Ueffect,input,B1,B2,bc] = InitWFSim(Wp,options,plotMesh);

if Animate > 0
    scrsz = get(0,'ScreenSize');
    hfig = figure('color',[0 166/255 214/255],'units','normalized','outerposition',...
           [0 0 1 1],'ToolBar','none','visible', 'on'); 
end

%% Loop 
for k=1:Wp.sim.NN
    tic
    it        = 0;
    eps       = 1e19;
    epss      = 1e20;

    while ( eps>conv_eps && it<max_it && eps<epss );
        it   = it+1;
        epss = eps;
        
        if k>1 
            max_it = max_it_dyn; 
        end
        
%         if k>=20
%             input{k}.beta = input{k}.beta+[0;.1];
%         end
        if k>=70
            input{k}.beta = input{k}.beta+[0;0];
        end


        
        [sys,Power(:,k),Ueffect(:,k),a(:,k),CT(:,k),Wp] = ...
                    Make_Ax_b(Wp,sys,sol,input{k},B1,B2,bc,k,options);              % Create system matrices
        [sol,sys] = Computesol(sys,input{k},sol,k,it,options);                      % Compute solution
        [sol,eps] = MapSolution(Wp.mesh.Nx,Wp.mesh.Ny,sol,k,it,options);            % Map solution to field
        
    end
 
    if Animate > 0
        if ~rem(k,Animate)
            Animation; 
        end; 
    end; 
    

end

%%
figure(2);clf
subplot(3,2,1)
plot(Power(1,:));grid;xlabel('k');ylabel('P_1')
subplot(3,2,3)
plot(a(1,:));grid;xlabel('k');ylabel('a_1')
subplot(3,2,5)
plot(Ueffect(1,:));grid;xlabel('k');ylabel('U_1')

subplot(3,2,2)
plot(Power(2,:));grid;xlabel('k');ylabel('P_2')
subplot(3,2,4)
plot(a(2,:));grid;xlabel('k');ylabel('a_2')
subplot(3,2,6)
plot(Ueffect(2,:));grid;xlabel('k');ylabel('U_2')


