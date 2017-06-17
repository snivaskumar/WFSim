function [ Wp, input ] = meshing( Wp, plotMesh, PrintGridMismatch )

% Default settings
if nargin <= 0; error('Please specify a meshing case.'); end;
if nargin <= 1; plotMesh = true;                         end;
if nargin <= 2; PrintGridMismatch = true;                end;

% Some pre-processing to determine correct input file location
meshingloc  = which('meshing.m');
strbslash   = strfind(meshingloc,'\');
WFSimfolder = meshingloc(1:strbslash(end-2));

switch lower(Wp.name)
    
    case lower('RobustMpc')
        type   = 'lin';          % Meshing type ('lin' or 'exp')
        Lx     = 2000;           % Domain length in x-direction (m)
        Ly     = 600;           % Domain length in y-direction (m)
        Nx     = 50;             % Number of grid points in x-direction
        Ny     = 25;             % Number of grid points in y-direction
        Crx    = [400+6*90];     % Turbine locations in x-direction (m)
        Cry    = [300];       % Turbine locations in y-direction (m)
        
        loadedinput = load([WFSimfolder 'Data_SOWFA\YawCase3\system_input.mat']); % load input settings
        loadedinput.input.phi = 0*loadedinput.input.phi;
        
        % Correctly format inputs (temporary function)
        for j = 1:length(loadedinput.input.t)
            input{j}.t    = loadedinput.input.t(j);
            input{j}.beta = [loadedinput.input.beta(j,1)'];
            input{j}.phi  = [loadedinput.input.phi(j,1)'];
        end;
        
        % Calculate delta inputs
        for j = 1:length(loadedinput.input.t)-1
            input{j}.dbeta = [loadedinput.input.beta(j+1,1)'- loadedinput.input.beta(j,1)'];
            input{j}.dphi  = [loadedinput.input.phi(j+1,1)' - loadedinput.input.phi(j,1)'] ;
        end;
        
        Drotor      = 90;     % Turbine rotor diameter in (m)
        powerscale  = 1.0;    % Turbine powerscaling
        forcescale  = 1;    % Turbine force scaling
        
        h        = 1.0;       % Sampling time (s)
        L        = 350;       % Simulation length (s)
        mu       = 0*18e-5;     % Dynamic flow viscosity
        Rho      = 1.20;      % Flow density (kg m-3)
        u_Inf    = 9.0;       % Freestream flow velocity x-direction (m/s)
        v_Inf    = 0.0;       % Freestream flow velocity y-direction (m/s)
        p_init   = 0.0;       % Initial values for pressure terms (Pa)
        
        lmu      = 1.75;         % Mixing length in x-direction (m)
        turbul   = true;      % Use mixing length turbulence model (true/false)
        n        = 2;
        m        = 6;
        
        
    otherwise
        error('No valid meshing specified. Please take a look at Wp.name.');
end;

%% Calculate time
time = 0:h:L;          % time vector for simulation
NN   = length(time)-1; % Total number of simulation steps

if strcmp(lower(type),'lin')
    % linear gridding
    ldx  = linspace(0,Lx,Nx);
    ldy  = linspace(0,Ly,Ny);
else
    error('Wrong meshing type specified in "type".');
end;

Nx   = length(ldx);
Ny   = length(ldy);
ldxx = repmat(ldx',1,Ny);
ldyy = repmat(ldy,Nx,1);

% Create secondary grid from primary grid
ldx2  = 0.5*(ldx(1:end-1)+ldx(2:end));
ldx2  = [ldx2 2*ldx2(end)-ldx2(end-1)]; % add extra cells
ldy2  = 0.5*(ldy(1:end-1)+ldy(2:end));
ldy2  = [ldy2 2*ldy2(end)-ldy2(end-1)]; % add extra cells
ldxx2 = repmat(ldx2',1,Ny);
ldyy2 = repmat(ldy2,Nx,1);

% Calculate cell dimensions
dx   = diff(ldx);
dxx  = repmat([dx'; dx(end)],1,Ny);
dx2  = diff(ldx2);
dxx2 = repmat([dx2'; dx2(end)],1,Ny);
dy   = diff(ldy);
dyy  = repmat([dy, dy(end)],Nx,1);
dy2  = diff(ldy2);
dyy2 = repmat([dy2, dy2(end)],Nx,1);

% Calculate location of turbines in grid and grid mismatch
for i = 1:length(Crx)
    % Calculate cells relevant for turbine (x-dir) on primary grid
    [~,xline(i,1)] = min(abs(ldx-Crx(i))); % automatically picks earliest entry in vector
    
    % Calculate cells closest to turbines (y-dir) on both grids
    [ML_prim, L_prim ] = min(abs(ldy- (Cry(i)-Drotor/2)));
    [ML_sec , L_sec  ] = min(abs(ldy2-(Cry(i)-Drotor/2)));
    [MR_prim, R_prim ] = min(abs(ldy- (Cry(i)+Drotor/2)));
    [MR_sec , R_sec  ] = min(abs(ldy2-(Cry(i)+Drotor/2)));
    
    yline{i}  = L_prim:1: R_prim; % turbine cells for primary grid
    % ylinev{i} = L_sec :1: R_sec ; % turbine cells for secondary grid
    ylinev{i} = L_prim:1: R_prim+1; % JWs code -- Bart: I feel like this needs fixing
    
    if PrintGridMismatch
        % Calculate turbine-grid mismatch
        disp([' TURBINE ' num2str(i) ' GRID MISMATCH:']);
        disp(['                    Primary           Secondary ']);
        disp(['       center:   (' num2str(min(abs(Crx(i)-ldx)),'%10.2f\n') ',' num2str(min(abs(Cry(i)-ldy)),'%10.2f\n') ') m.      (' num2str(min(abs(Crx(i)-ldx2)),'%10.2f\n') ',' num2str(min(abs(Cry(i)-ldy2)),'%10.2f\n') ') m.']);
        disp(['   left blade:   (' num2str(min(abs(Crx(i)-ldx)),'%10.2f\n') ',' num2str(ML_prim,'%10.2f\n')             ') m.      (' num2str(min(abs(Crx(i)-ldx2)),'%10.2f\n') ',' num2str(ML_sec,'%10.2f\n')                ') m.']);
        disp(['  right blade:   (' num2str(min(abs(Crx(i)-ldx)),'%10.2f\n') ',' num2str(MR_prim,'%10.2f\n')             ') m.      (' num2str(min(abs(Crx(i)-ldx2)),'%10.2f\n') ',' num2str(MR_sec,'%10.2f\n')                ') m.']);
        disp(' ');
    end;
end;

% Calculate state sizes
Nu = (Nx-3)*(Ny-2);   % Number of u velocities in state vector
Nv = (Nx-2)*(Ny-3);   % Number of v velocities in state vector
Np = (Nx-2)*(Ny-2)-2; % Number of pressure terms in state vector


%% Export to Wp and input
Wp         = struct('Nu',Nu,'Nv',Nv,'Np',Np,'name',Wp.name);
Wp.sim     = struct('h',h,'time',time,'L',L,'NN',NN);
Wp.turbine = struct('Drotor',Drotor,'powerscale',powerscale,'forcescale',forcescale, ...
    'N',length(Crx),'Crx',Crx,'Cry',Cry);
Wp.site    = struct('mu',mu,'Rho',Rho,'u_Inf',u_Inf,'v_Inf',v_Inf,'p_init',p_init, ...
    'lmu',lmu,'turbul',turbul,'m',m,'n',n);
Wp.mesh    = struct('Lx',Lx,'Ly',Ly,'Nx',Nx,'Ny',Ny,'ldxx',ldxx,'ldyy',ldyy,'ldxx2',...
    ldxx2,'ldyy2',ldyy2,'dxx',dxx,'dyy',dyy,'dxx2',dxx2,'dyy2',dyy2,...
    'xline',xline,'type',type);Wp.mesh.yline = yline; Wp.mesh.ylinev = ylinev; % Do not support struct command
end