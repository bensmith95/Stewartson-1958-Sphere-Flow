%% Solves the Stewartson equations for the Impinging Region 
% Inputs:
% Re - Reynolds Number (sqrt) 

function Stewartson(Re)
    %% Initialise
    str1 = fprintf('Initialising...\n');
    %  set step size
    h = 2^-4;
    % set eta
    eta = 0:h:30; Neta = length(eta); 
    % set beta
    beta = 0:h:10; Nbeta = length(beta); 

    % assign variable storage
    Psi = zeros(Nbeta+2,Neta+2); Omega = zeros(Nbeta+2,Neta+2); 

    % Obtain and update BC's
    hi = 1/h; [Pxxt,Pt] = BC(h,Re);
    % Psi
    Psi(Nbeta+1,2:Neta+1) = Pt(:)'; 
    % W
    Omega(Nbeta+1,2:Neta+1) = 2*(Psi(Nbeta+1,2:Neta+1)-Psi(Nbeta,2:Neta+1))*hi^2 - Pxxt(:)';
    Omega(2:Nbeta+1,2) = -2*Psi(2:Nbeta+1,3)*hi^2;
    Omega(3:Nbeta,Neta+1) = (2*Psi(3:Nbeta,Neta+1)-Psi(2:Nbeta-1,Neta)-Psi(4:Nbeta+1,Neta))*hi^2;
    % Update extrapolated points
    Psi(2:Nbeta+1,1) = Psi(2:Nbeta+1,3);
    Psi(2:Nbeta+1,Neta+2) = Psi(2:Nbeta+1,Neta);
    Psi(1,2:Neta+1) = 3*Psi(2,2:Neta+1) - 3*Psi(3,2:Neta+1) + Psi(4,2:Neta+1);
    Psi(Nbeta+2,2:Neta+1) = Psi(Nbeta,2:Neta+1);
    Omega(2:Nbeta+1,1) = 3*Omega(2:Nbeta+1,2) - 3*Omega(2:Nbeta+1,3) + Omega(2:Nbeta+1,4);
    Omega(2:Nbeta+1,Neta+2) = 3*Omega(2:Nbeta+1,Neta+1) - 3*Omega(2:Nbeta+1,Neta) + Omega(2:Nbeta+1,Neta-1);
    Omega(1,2:Neta+1) = -3*Omega(3,2:Neta+1) + Omega(4,2:Neta+1);
    Omega(Nbeta+2,2:Neta+1) = 3*Omega(Nbeta+1,2:Neta+1) - 3*Omega(Nbeta,2:Neta+1) + Omega(Nbeta-1,2:Neta+1);

    %% FMG-FAS Alg.
    fprintf(repmat('\b',1,str1)); str1 = fprintf('FMG-FAS Algorithm...\n');
    [Psi,Omega,~] = FMG(Psi,Omega,@BC,h,Re); 
    fprintf(repmat('\b',1,str1)); str1 = fprintf('Solution Converged.\n');

    %% POST PROCESS
    str2 = fprintf('Post-Processing...\n');
    % load Boundary Layer flow
    filename = 'BL.mat';
    load(filename,'VelBL','eta','theta'); UB = VelBL{1}; Eta = 0:h:30; 
    % Determine U inlet
    Tm = -10/Re + pi/2; [~,i] = min(abs(theta-Tm));
    UI = spline(eta,-UB(i,:),Eta); % U in
    eta = Eta;
    % Determine W via finite differences of SF
    W = zeros(Nbeta,Neta); 
    W(2:Nbeta-1,:) = (Psi(4:Nbeta+1,2:Neta+1)-Psi(2:Nbeta-1,2:Neta+1))*0.5*hi;
    W(1,:) = -(3*Psi(2,2:Neta+1)-4*Psi(3,2:Neta+1)+Psi(4,2:Neta+1))*0.5*hi;
    W(Nbeta,:) = 0; W(:,1) = 0;  
    % Determine U via finite differences of SF
    U = zeros(Nbeta,Neta);
    U(:,2:Neta-1) = -(Psi(2:Nbeta+1,4:Neta+1)-Psi(2:Nbeta+1,2:Neta-1))*0.5*hi;
    U(:,1) = 0; U(:,Neta) = 0; U(Nbeta,:) = UI; U(1,:) = 0;
    
    % Determine V
    [V,~] = Vvel(W,U,Re,h); 
    
    % Solve Poisson Equation for Pressure
    etam = (eta(1:end-1)+eta(2:end))*0.5; betam = (beta(1:end-1)+beta(2:end))*0.5;
    P = Pressure(W,U,Re,h); 

    % save data to file
    VelIMP{1} = U; VelIMP{2} = V; VelIMP{3} = W;
    VelIMP{4} = Psi; VelIMP{5} = Omega; VelIMP{6} = P;
    filename = ['Stew_Re=',num2str(Re),'.mat'];
    save(filename, 'VelIMP', 'eta', 'beta', 'etam', 'betam')
    fprintf(repmat('\b',1,str2)); str2 = fprintf('Flow saved in %s\n', filename); pause(1)

    fprintf(repmat('\b',1,str2)); fprintf(repmat('\b',1,str1));
end
%% BOUNDARY CONDITIONS
function [Pxx,Psit] = BC(h,Re)
    hi = 1/h;
    % load Boundary Layer flow 
    filename = 'BL.mat';
    load(filename,'VelBL','eta','theta'); UB = VelBL{1}; Eta = 0:h:30; 
    % Determine betamax in theta 
    Tm = -10/Re + pi/2; [~,i] = min(abs(theta-Tm));
    % Determine U inlet
    UI = spline(eta,-UB(i,:),Eta); % U in
    % Vorticity BC's
    Pxx = zeros(size(Eta)); Pxx(2:end-1) = (UI(3:end)-UI(1:end-2))*0.5*hi;
    Pxx(1) = -(3*UI(1)-4*UI(2)+UI(3))*0.5*hi; Pxx(end) = (3*UI(end)-4*UI(end-1)+UI(end-2))*0.5*hi; Pxx=-Pxx;
    % SF BC's
    Psit = -cumtrapz(Eta,UI); % Psi Top
end
