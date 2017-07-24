clear;clc;close all;

WFSim_addpaths

%% Initialize script
options.Projection     = 0;                      % Use projection (true/false)
options.Linearversion  = 0;                      % Provide linear variant of WFSim (true/false)
options.exportLinearSol= 0;                      % Calculate linear solution of WFSim
options.Derivatives    = 0;                      % Compute derivatives
options.startUniform   = 0;                      % Start from a uniform flowfield (true) or a steady-state solution (false)
options.exportPressures= ~options.Projection;    % Calculate pressure fields

Wp.name             = 'RobustMpc';      % Meshing name (see "\bin\core\meshing.m")
Wp.Turbulencemodel  = 'WFSim3';

Animate       = 1;                      % Show 2D flow fields every x iterations (0: no plots)
plotMesh      = 0;                      % Show meshing and turbine locations
conv_eps      = 1e-6;                   % Convergence threshold
max_it_dyn    = 1;                      % Maximum number of iterations for k > 1

if options.startUniform==1;max_it = 1;else max_it = 50;end

[Wp,sol,sys,Power,CT,a,Ueffect,input,B1,B2,bc] = InitWFSim(Wp,options,plotMesh);
if Animate > 0
    scrsz = get(0,'ScreenSize');
    hfig = figure('color',[0 166/255 214/255],'units','normalized','outerposition',...
        [0 0 1 1],'ToolBar','none','visible', 'on');
end
% Loop converges steady-stte
for k=1:1
    tic
    it        = 0;
    eps       = 1e19;
    epss      = 1e20;
    while ( eps>conv_eps && it<max_it && eps<epss );
        it   = it+1;
        epss = eps;
        if k>1;max_it = max_it_dyn;end
        [sys,Power(:,k),Ueffect(:,k),a(:,k),CT(:,k)] = ...
            Make_Ax_b(Wp,sys,sol,input{k},B1,B2,bc,k,options);              % Create system matrices
        [sol,sys] = Computesol(sys,input{k},sol,k,it,options);                      % Compute solution
        [sol,eps] = MapSolution(Wp.mesh.Nx,Wp.mesh.Ny,sol,k,it,options);            % Map solution to field
    end
    if Animate>0;if~rem(k,Animate);Animation;end;end;
end;
disp('Wind farm in steady-state.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize and apply perturbations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ny1           = Wp.mesh.yline{1};
nx1           = Wp.mesh.xline(1);
ny2           = Wp.mesh.yline{2};
nx2           = Wp.mesh.xline(2);
ny3           = Wp.mesh.yline{3};
nx3           = Wp.mesh.xline(3);
uss           = sol.u; % steady-state solution

Q             = 10; % number of perturbations
perturbation  = 2*rand(Q,length(Wp.mesh.yline{1}))-1;

% Time loop
for l=1:Q
    
    % Apply perturbation in front of first turbine
    sol.u(nx1-1,ny1) = uss(nx1-1,ny1) + 3*perturbation(l,:);
    
    for k=2:Wp.sim.NN
        
        % Save velocities in front and behind turbines
        flow.T1.up(:,:,k,l)   = sol.u(nx1-1,ny1);  % flow.turbine1.upwind(x,y,time,perturbation)
        flow.T1.down(:,:,k,l) = sol.u(nx1+1,ny1);  % flow.turbine2.downwind(x,y,time,perturbation)
        flow.T2.up(:,:,k,l)   = sol.u(nx2-1,ny2);  % maybe take mean here
        flow.T2.down(:,:,k,l) = sol.u(nx2+1,ny2);
        flow.T3.up(:,:,k,l)   = sol.u(nx3-1,ny3);
        flow.T3.down(:,:,k,l) = sol.u(nx3+1,ny3);
        
        % Build and iterate wind farm
        [sys,Power(:,k),Ueffect(:,k),a(:,k),CT(:,k),Wp] = ...
            Make_Ax_b(Wp,sys,sol,input{k},B1,B2,bc,k,options);              % Create system matrices
        [sol,sys] = Computesol(sys,input{k},sol,k,it,options);              % Compute solution
        [sol,eps] = MapSolution(Wp.mesh.Nx,Wp.mesh.Ny,sol,k,it,options);    % Map solution to field
                
        %Plot
        if Animate>0;if~rem(k,Animate);Animation;end;end;
        
    end
end
%%
figure(2);clf
subplot(3,3,1)
plot(Power(1,:));grid;xlabel('k');ylabel('P_1')
subplot(3,3,4)
plot(a(1,:));grid;xlabel('k');ylabel('a_1')
subplot(3,3,7)
plot(Ueffect(1,:));grid;xlabel('k');ylabel('U_1')

subplot(3,3,2)
plot(Power(2,:));grid;xlabel('k');ylabel('P_2')
subplot(3,3,5)
plot(a(2,:));grid;xlabel('k');ylabel('a_2')
subplot(3,3,8)
plot(Ueffect(2,:));grid;xlabel('k');ylabel('U_2')

subplot(3,3,3)
plot(Power(3,:));grid;xlabel('k');ylabel('P_3')
subplot(3,3,6)
plot(a(3,:));grid;xlabel('k');ylabel('a_3')
subplot(3,3,9)
plot(Ueffect(3,:));grid;xlabel('k');ylabel('U_3')
