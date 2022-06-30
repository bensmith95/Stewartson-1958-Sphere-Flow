%% Obtain V velocity component
% Inputs:
% W - W velocity component      U - U velocity component
% h - size of grid              Re - Reynolds Number
% Outputs:
% V - V velocity component      r - final residual
function [V,r] = Vvel(W,U,Re,h)
    %% Initalise
    str1 = fprintf('Calculating V...\n');
    [Ny,Nx] = size(U); V = zeros(Ny+2,Nx+2); [Vt,Vl] = VBCs(h,Re);
    % set Dirichlet Bcs
    V(2:Ny+1,2) = Vl;  
    V(Ny+1,2:Nx+1) = Vt; 
    % set Extrapolated Points
    V(2:Ny+1,1) = 3*V(2:Ny+1,2) - 3*V(2:Ny+1,3) + V(2:Ny+1,4);
    V(Ny+2,2:Nx+1) = 3*V(Ny+1,2:Nx+1) - 3*V(Ny,2:Nx+1) + V(Ny-1,2:Nx+1);

    %% V-cycles until convergence
    F = zeros(size(V)); tol = 1e-5; r = Inf;
    str2 = fprintf('residual = %f, %.2f%% Complete\n',[NaN,0]);
    while r>tol 
        fprintf(repmat('\b',1,str2))
        str2 = fprintf('residual = %f, %.2f%% Complete\n',[r,tol/r*100]);
        [V,r] = VCycle(V,W,U,F,Re,h,1); 
    end
    fprintf(repmat('\b',1,str2)); fprintf(repmat('\b',1,str1));

    % Truncate V to size of domain
    V = V(2:Ny+1,2:Nx+1);
end
%% DEFINE BOUNDARY CONDITIONS
function [Vt,Vl] = VBCs(h,Re)
    % load Boundary Layer flow 
    filename = 'BL.mat';
    load(filename,'VelBL','eta','theta'); VB = VelBL{2}; Eta = eta; eta = 0:h:30; 
    % Determine betamax in theta
    Tm = -10/Re + pi/2; [~,i] = min(abs(theta-Tm));
    % map to current grid
    VI = spline(Eta,VB(i,:),eta); % V in
    % set BCs
    Vt = VI; Vl = 1;
end
%% V-CYCLE 
function [V,r] = VCycle(V,W,U,F,Re,h,w)
	% Pre-Smoothing
	V = SOR(V,W,U,F,Re,h,5,w);
	% Compute Residual
	[res,~] = Vresidual(V,W,U,F,Re,h);
	% Restrict to coarse grid
	rhs = Vrestrict(res); W2h = W(1:2:end,1:2:end); U2h = U(1:2:end,1:2:end); 
    % guess for error
	err = zeros(size(rhs));
	% stop recursion at smallest grid size, otherwise continue recursion
    if h>0.8
        	err = SOR(err,W2h,U2h,rhs,Re,2*h,10,1); 
	else        
        	[err,~] = VCycle(err,W2h,U2h,rhs,Re,2*h,w);    
    end
	% Prolongation and Correction
	V = V + prolong(err); N = size(V); Nx = N(2)-2; Ny = N(1)-2;
    % set Extrapolated Points
    V(2:Ny+1,1) = 3*V(2:Ny+1,2) - 3*V(2:Ny+1,3) + V(2:Ny+1,4);
    V(Ny+2,2:Nx+1) = 3*V(Ny+1,2:Nx+1) - 3*V(Ny,2:Nx+1) + V(Ny-1,2:Nx+1);
    % set neumann points
    V(1,3:Nx+1) = V(3,3:Nx+1);
    V(2:Ny,Nx+2) = V(2:Ny,Nx);
	% Post-Smoothing
	V = SOR(V,W,U,F,Re,h,3,w); 
    % max residual of correction
    [~,r] = Vresidual(V,W,U,F,Re,h);
end
%% Succesive Over-Relaxation
function V = SOR(V,W,U,S,Re,h,n,w)
    [Ny,Nx] = size(U); hi = 1/h; eps = 1/Re; 
    for k=1:n
        for j=Ny:-1:2
            for i=3:Nx+1
                % W Upwind
                if W(j-1,i-1)<=0
                    if i==Nx+1; a=0; CW=0; else; a = -3*0.5*hi; CW = -(-4*V(j,i+1)+V(j,i+2))*0.5*hi; end
                else
                    if i==Nx+1; a=0; CW=0; else; a = 3*0.5*hi; CW = (-4*V(j,i-1)+V(j,i-2))*0.5*hi; end
                end
                % U Upwind
                if U(j-1,i-1)<=0
                    b = -3*0.5*hi; CU = -(-4*V(j+1,i)+V(j+2,i))*0.5*hi;
                else
                    b = 3*0.5*hi; CU = (-4*V(j-1,i)+V(j-2,i))*0.5*hi;
                end
                % Solve Linear Equation
                N = S(j,i) + eps*(V(j-1,i)+V(j,i-1)+V(j,i+1)+V(j+1,i))*hi^2 - W(j-1,i-1)*CW - U(j-1,i-1)*CU;
                D = W(j-1,i-1)*a + U(j-1,i-1)*b + 4*eps*hi^2;
                VT = N/D;
                % Under-Relax
                V(j,i) = (1-w)*V(j,i) + w*VT;
                % Update Extrap. points
                if ismember(i,[3,4]); V(j,1) = 3*V(j,2) - 3*V(j,3) + V(j,4); end
                if ismember(j,[Ny-1,Ny]); V(Ny+2,i) = 3*V(Ny+1,i) - 3*V(Ny,i) + V(Ny-1,i); end
                % Update Neumann BC points
                if j==3 && ismember(i,3:Nx+1); V(1,i) = V(3,i); end 
                if i==Nx && ismember(j,2:Ny); V(j,Nx+2) = V(j,Nx); end
            end
        end
    end
end
%% Residual
function [res,r] = Vresidual(V,W,U,S,Re,h)
    [Ny,Nx] = size(U); hi = 1/h; eps = 1/Re; res = zeros(Ny+2,Nx+2);
    for j=2:Ny
        for i=3:Nx
            if W(j-1,i-1)<=0
                if i==Nx+1; dVdeta = 0; else; dVdeta = -(3*V(j,i)-4*V(j,i+1)+V(j,i+2))*0.5*hi; end
            else
                if i==Nx+1; dVdeta = 0; else; dVdeta = (3*V(j,i)-4*V(j,i-1)+V(j,i-2))*0.5*hi; end
            end
            if U(j-1,i-1)<=0
                dVdbeta = -(3*V(j,i)-4*V(j+1,i)+V(j+2,i))*0.5*hi;
            else
                dVdbeta = (3*V(j,i)-4*V(j-1,i)+V(j-2,i))*0.5*hi;
            end
            LapV = (V(j-1,i)+V(j,i-1)-4*V(j,i)+V(j,i+1)+V(j+1,i))*hi^2;
            res(j,i) = S(j,i) + eps*LapV - W(j-1,i-1)*dVdeta - U(j-1,i-1)*dVdbeta;
        end
    end
    r = max(abs(res),[],'all');
end
%% Restrict Residual to Coarse grid
function r2h = Vrestrict(rh)
    N = 0.5*(size(rh)-3)+1; r2h = zeros(N+2);

    for j = 3:N(1)
        for i = 3:N(2)
            r2h(j,i) = 0.0625*(rh(2*j-3,2*i-3)+rh(2*j-3,2*i-1)+rh(2*j-1,2*i-3)+rh(2*j-1,2*i-1)) + ...
                     0.125*(rh(2*j-2,2*i-3)+rh(2*j-2,2*i-1)+rh(2*j-3,2*i-2)+rh(2*j-1,2*i-2)) + ...
                     0.25*rh(2*j-2,2*i-2);
        end
    end
    r2h([2,N(1)+1],2:N(2)+1) = rh([2,end-1],2:2:end-1); r2h(2:N(1)+1,[2,N(2)+1]) = rh(2:2:end-1,[2,end-1]);
end
%% Interpolate to Fine grid
function F = prolong(C)
    N = size(C); F = zeros(2*(N-3)+1); M = size(F)+2; F = zeros(M);

    F(2:2:M(1)-1,2:2:M(2)-1) = C(2:N(1)-1,2:N(2)-1);
    F(2:2:M(1)-1,3:2:M(2)-2) = 0.5*(C(2:N(1)-1,2:N(2)-2)+C(2:N(1)-1,3:N(2)-1));
    F(3:2:M(1)-2,2:2:M(2)-1) = 0.5*(C(2:N(1)-2,2:N(2)-1)+C(3:N(1)-1,2:N(2)-1));  
    F(3:2:M(1)-2,3:2:M(2)-2) = 0.25*(C(2:N(1)-2,2:N(2)-2)+C(2:N(1)-2,3:N(2)-1) + ...
                                     C(3:N(1)-1,2:N(2)-2)+C(3:N(1)-1,3:N(2)-1));
end
