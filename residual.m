%% Residual of the SF-Vort. equations 
% i.e. how well the solution satisfies the system
% Inputs:
% Psi - Approximation for SF    W - Approximation for Vort.
% Sp - RHS of SF equation       Sw - RHS of Vort. eqn
% h - size of grid              Re - Reynolds Number          
% Outputs:
% Pres - residual of SF eqn     Wres - residual of Vort. eqn 
function [Pres,Wres] = residual(Psi,W,Sp,Sw,h,Re)
    %% Initialise
    N = size(Psi); Nx = N(2)-2; Ny = N(1)-2;
    Pres = zeros(N); Wres = zeros(N);
    eps = 1/Re; hi = 1/h;
    
    %% Iterate
    for j = Ny:-1:3
        for i = 3:Nx+1
            dPsidx = (Psi(j,i+1)-Psi(j,i-1))*0.5*hi;
            dPsidy = (Psi(j+1,i)-Psi(j-1,i))*0.5*hi;
            if i<Nx+1
                % SF
                if dPsidy>=0
                    dWdx = (3*W(j,i)-4*W(j,i-1)+W(j,i-2))*0.5*hi; 
                else
                    dWdx = -(3*W(j,i)-4*W(j,i+1)+W(j,i+2))*0.5*hi; 
                end
                if -dPsidx>=0
                    dWdy = (3*W(j,i)-4*W(j-1,i)+W(j-2,i))*0.5*hi; 
                else
                    dWdy = -(3*W(j,i)-4*W(j+1,i)+W(j+2,i))*0.5*hi; 
                end
                LapW = (W(j-1,i)+W(j,i-1)-4*W(j,i)+W(j,i+1)+W(j+1,i))*hi^2;
                Pres(j,i) = Sp(j,i) - (dPsidy*dWdx - dPsidx*dWdy - eps*LapW);

                % Vort.
                LapPsi = (Psi(j-1,i)+Psi(j,i-1)-4*Psi(j,i)+Psi(j,i+1)+Psi(j+1,i))*hi^2;
                Wres(j,i) = Sw(j,i) - (W(j,i) + LapPsi);
            else
                Pres(j,Nx+1) = Sp(j,Nx+1) - (Psi(j-1,Nx+1)-Psi(j-1,Nx)+Psi(j,Nx)-2*Psi(j,Nx+1)+Psi(j,Nx+2)+Psi(j+1,Nx+1)-Psi(j+1,Nx));
            end
        end
    end
end