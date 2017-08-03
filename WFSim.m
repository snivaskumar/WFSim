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

Animate       = 0;                      % Show 2D flow fields every x iterations (0: no plots)
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
clear Power a Ueffect Wp.sim.time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize and apply perturbations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ny1           = Wp.mesh.yline{1};
nx1           = Wp.mesh.xline(1);
ny2           = Wp.mesh.yline{2};
nx2           = Wp.mesh.xline(2);
ny3           = Wp.mesh.yline{3};
nx3           = Wp.mesh.xline(3);
uss           = sol.u; % steady-state solutions
vss           = sol.v;
pss           = sol.p;
uInfss        = sol.u(1); % steady-state free-stream velocity

% Number of perturbed trajectories
Ns            = 5;
% Number of samples in one time series
Np            = 100;
% Vf is one measurement trajectory
Vf            = repmat(uInfss+2-4*rand(1,Np),Ns,1);%repmat(uInfss+.5*[zeros(1,100) ones(1,400)],Ns,1);%
% Vs is random perturbed measurement trajectory
Vs            = Vf + 2*rand(Ns,Np);
% Time vector simulation
Wp.sim.h      = 1;  
Wp.sim.time   = 0:Wp.sim.h:Np;

for l=1:Ns
    
    % Set to steady-state
    sol.u                   = uss;
    sol.v                   = vss;
    sol.p                   = pss;
    
    
    % Time loop
    for k=1:Np
        
        % Apply perturbation
        Wp.site.u_Inf           = Vs(l,k);
        [B1,B2,bc]              = Compute_B1_B2_bc(Wp);
        sol.u(1:2,1:Wp.mesh.Ny) = Wp.site.u_Inf;
        
        
        % Save velocities in front and behind turbines
        flow.T0.up(k,l)   = mean(sol.u(1,ny1));      % flow upwind windfarm
        flow.T1.up(k,l)   = mean(sol.u(nx1-1,ny1));  % flow.turbine1.upwind(x,y,time,perturbation)
        flow.T1.down(k,l) = mean(sol.u(nx1+1,ny1));  % flow.turbine2.downwind(x,y,time,perturbation)
        flow.T2.up(k,l)   = mean(sol.u(nx2-1,ny2));
        flow.T2.down(k,l) = mean(sol.u(nx2+1,ny2));
        flow.T3.up(k,l)   = mean(sol.u(nx3-1,ny3));
        flow.T3.down(k,l) = mean(sol.u(nx3+1,ny3));
        
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
% Turbine transfer from w = [dU da]^T to y = [dF dP]^T
y1 = Wp.turbine.yT{1};
y2 = Wp.turbine.yT{2};
y3 = Wp.turbine.yT{3};
w1 = Wp.turbine.wT{1};
w2 = Wp.turbine.wT{2};
w3 = Wp.turbine.wT{3};

figure(2);clf
subplot(3,3,1)
plot(y1(2,:));grid;xlabel('k');ylabel('\delta P_1');axis tight
subplot(3,3,4)
plot(w1(2,:));grid;xlabel('k');ylabel('\delta a_1');axis tight
subplot(3,3,7)
plot(w1(1,:));grid;xlabel('k');ylabel('\delta U_1');axis tight
subplot(3,3,2)
plot(y2(2,:));grid;xlabel('k');ylabel('\delta P_2');axis tight
subplot(3,3,5)
plot(w2(2,:));grid;xlabel('k');ylabel('\delta a_2');axis tight
subplot(3,3,8)
plot(w2(1,:));grid;xlabel('k');ylabel('\delta U_2');axis tight
subplot(3,3,3)
plot(y3(2,:));grid;xlabel('k');ylabel('\delta P_3');axis tight
subplot(3,3,6)
plot(w3(2,:));grid;xlabel('k');ylabel('\delta a_3');axis tight
subplot(3,3,9)
plot(w3(1,:));grid;xlabel('k');ylabel('\delta U_3');axis tight

figure(3);clf
subplot(4,3,[1 3])
plot(flow.T0.up(:,l));grid;xlabel('k');ylabel('v_0');axis tight
subplot(4,3,4)
plot(flow.T1.up(:,l));grid;xlabel('k');ylabel('v_1^{in}');axis tight
subplot(4,3,7)
plot(flow.T1.down(:,l));grid;xlabel('k');ylabel('v_1^{out}');axis tight
subplot(4,3,5)
plot(flow.T2.up(:,l));grid;xlabel('k');ylabel('v_2^{in}');axis tight
subplot(4,3,8)
plot(flow.T2.down(:,l));grid;xlabel('k');ylabel('v_2^{out}');axis tight
subplot(4,3,6)
plot(flow.T3.up(:,l));grid;xlabel('k');ylabel('v_3^{in}');axis tight
subplot(4,3,9)
plot(flow.T3.down(:,l));grid;xlabel('k');ylabel('v_3^{out}');axis tight

