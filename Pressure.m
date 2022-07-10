%% Determine Pressure distribution
% via the Pressure Laplace equation
function P = Pressure(U,W,Re,h)
    %% Initialise
    [Nbeta,Neta] = size(U); hi = 1/h;
    
    % load Pressure inlet
    filename = 'BL.mat';
    load(filename,'VelBL','eta','theta'); PB = VelBL{4}; Eta = eta; eta = 0:h:30;
    % Determine betamax in theta
    Tm = -10/Re + pi/2; [~,i] = min(abs(theta-Tm));
    % map to current grid
    PI = spline(Eta,PB(i,:),eta); % P in
    
   % determine derivatives
    dUdeta = zeros(Nbeta,Neta); dUdbeta = dUdeta; 
    dWdeta = dUdeta; dWdbeta = dUdeta; 
    % U
    dUdeta(:,1) = -(3*U(:,1)-4*U(:,2)+U(:,3))*0.5*hi;
    dUdeta(:,2:Neta-1) = (U(:,3:Neta)-U(:,1:Neta-2))*0.5*hi; 
    dUdeta(:,Neta) = (3*U(:,Neta)-4*U(:,Neta-1)+U(:,Neta-2))*0.5*hi;
    dUdbeta(1,:) = -(3*U(1,:)-4*U(2,:)+U(3,:))*0.5*hi;
    dUdbeta(2:Nbeta-1,:) = (U(3:Nbeta,:)-U(1:Nbeta-2,:))*0.5*hi;
    dUdbeta(Nbeta,:) = (3*U(Nbeta,:)-4*U(Nbeta-1,:)+U(Nbeta-2,:))*0.5*hi;
    % W
    dWdeta(:,1) = -(3*W(:,1)-4*W(:,2)+W(:,3))*0.5*hi;
    dWdeta(:,2:Neta-1) = (W(:,3:Neta)-W(:,1:Neta-2))*0.5*hi; 
    dWdeta(:,Neta) = (3*W(:,Neta)-4*W(:,Neta-1)+W(:,Neta-2))*0.5*hi;
    dWdbeta(1,:) = -(3*W(1,:)-4*W(2,:)+W(3,:))*0.5*hi;
    dWdbeta(2:Nbeta-1,:) = (W(3:Nbeta,:)-W(1:Nbeta-2,:))*0.5*hi;
    dWdbeta(Nbeta,:) = (3*W(Nbeta,:)-4*W(Nbeta-1,:)+W(Nbeta-2,:))*0.5*hi;

    % determine f
    F = 2*Re*(dWdeta.*dUdbeta - dUdeta.*dWdbeta)*h^2; 
    F(Nbeta,:) = PI;
    

    %% Laplace operator 
    e = ones(Neta,1); 
    D = spdiags([e,-4*e,e],-1:1,Neta,Neta); D(1,2) = 2; D(Neta,Neta-1) = 2;
    I = speye(Neta);
    Z = sparse(Neta,Neta); 
    A = [D, 2*I, repmat(Z,[1,Nbeta-2])];
    for j = 2:Nbeta-1
        A = [A; repmat(Z,[1,j-2]), I, D, I, repmat(Z,[1,Nbeta-3-j+2])];
    end
    A = [A; repmat(Z,[1,Nbeta-1]), I];

    % reorganise B into a vector
    f = zeros(Neta*Nbeta,1);
    for j=1:Nbeta
       f(1+(j-1)*Neta:j*Neta) = F(j,:); 
    end

    %% Solve linear system
    Pv = A\f;

    % reorganise solution into matrix form
    P = zeros(Nbeta,Neta);
    for j=1:Nbeta
        P(j,:) = Pv(1+(j-1)*Neta:j*Neta); 
    end
end
