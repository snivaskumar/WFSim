% Define spatial varying mixing-length parameter
xline            = sort(unique(xline));

if N==1
    x                = [zeros(1,xline(1)+n) linspace(0,lmu,Nx-xline(1)-n)]';
    y                = [zeros(1,yline{1}(1)-1) ones(1,length(yline{1})) zeros(1,Ny-yline{1}(end))] ;
    mixing_length    = (repmat(x,1,Ny).*repmat(y,Nx,1))*0.5*Drotor;
elseif N==2
    x                = [zeros(1,xline(1)+n) linspace(0,lmu,xline(1+1)-xline(1)-4*n)]';
    x                = [x;zeros(4*n,1);linspace(0,lmu,Nx-xline(2)-n)'];
    y                = [zeros(1,yline{1}(1)-1) ones(1,length(yline{1})) zeros(1,Ny-yline{1}(end))] ;
    mixing_length    = (repmat(x,1,Ny).*repmat(y,Nx,1))*0.5*Drotor;
    
elseif N==3
    xline            = sort(unique(xline));
    x                = [zeros(1,xline(1)+n) linspace(0,lmu,xline(2)-xline(1)-m)]';
    x                = [x;zeros(m,1);linspace(0,lmu,xline(3)-xline(2)-m)'];
    x                = [x;zeros(m,1);linspace(0,lmu,Nx-xline(3)-n)'];
    y                = [zeros(1,yline{1}(1)-1) ones(1,length(yline{1})) zeros(1,Ny-yline{1}(end))] ;
    mixing_length    = (repmat(x,1,Ny).*repmat(y,Nx,1))*0.5*Drotor;
end
H                = fspecial('disk',.5); % You need Nx,Nx to keep it symmetric
mixing_length    = filter2(H,mixing_length);

% Check with the following lines the mixing length parameter
%mixing_length(mixing_length<=10^-8)=nan;
%mesh(mixing_length);
%xlim([0 size(mixing_length,2)]);ylim([0 size(mixing_length,1)]);zlim([0 max(max(mixing_length))+.25*Drotor]);

% For u-momentum equation
ax.Tnx              = zeros(Nx,Ny);
ax.Tsx              = zeros(Nx,Ny);

ax.Tnx(2:Nx,1:Ny-1) = Rho*(mixing_length(2:Nx,1:Ny-1).^2).*(dxx(2:Nx,1:Ny-1)./(dyy(2:Nx,2:Ny).^2)).*abs(u(2:Nx,2:Ny)-u(2:Nx,1:Ny-1));
ax.Tsx(1:Nx-1,2:Ny) = Rho*(mixing_length(1:Nx-1,2:Ny).^2).*(dxx(2:Nx,2:Ny)./(dyy(2:Nx,2:Ny).^2)).*abs(u(2:Nx,1:Ny-1)-u(2:Nx,2:Ny));

ax.aN             = ax.aN + ax.Tnx;
ax.aS             = ax.aS + ax.Tsx;
ax.aP             = ax.aP + ax.Tnx + ax.Tsx;

% For v-momentum equation
ay.Tey            = zeros(Nx,Ny);
ay.Twy            = zeros(Nx,Ny);

ay.Tey(1:Nx-1,1:Ny) = Rho*(mixing_length(1:Nx-1,1:Ny).^2).*(dyy(1:Nx-1,1:Ny)./(dxx(1:Nx-1,1:Ny).^2)).*abs(v(2:Nx,1:Ny)-v(1:Nx-1,1:Ny));
ay.Twy(2:Nx,1:Ny)   = Rho*(mixing_length(2:Nx,1:Ny).^2).*(dyy(2:Nx,1:Ny)./(dxx(2:Nx,1:Ny).^2)).*abs(v(1:Nx-1,1:Ny)-v(2:Nx,1:Ny));

ay.aE             = ay.aE + ay.Tey;
ay.aW             = ay.aW + ay.Twy;
ay.aP             = ay.aP + ay.Tey + ay.Twy;