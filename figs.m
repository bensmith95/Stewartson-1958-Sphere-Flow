%% Impinging Region plots

clear all; clc; close all

% display info
set(0,'units','pixels'); disp = get(0,'ScreenSize');
Lx = disp(3); Ly = disp(4); 

%% Load flow
Re = 100;
IMPfile = ['Stew_Re=',num2str(Re),'.mat'];
load(IMPfile); U = VelIMP{1}; V = VelIMP{2}; W = VelIMP{3}; Psi = VelIMP{4}; Omega = VelIMP{5}; P = VelIMP{6};

%% Plot U-W Velocity Magnitude
figure(1); TT = ['Velocity Field: $\sqrt{R_e} =$ ',num2str(Re)];
fn = pcolor(eta,beta,sqrt(W.^2+U.^2)); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none'); 
fc = colorbar('location','southoutside'); caxis([0,0.15]); colormap('jet');
set(get(fc,'label'),'string','$\sqrt{W^2+U^2}$','interpreter','latex');
hold on; ii = 32; jj = 16; 
quiver(eta(1:ii:end),beta(1:jj:end),W(1:jj:end,1:ii:end),U(1:jj:end,1:ii:end),'color','k'); 
xlabel('\eta'); ylabel('-\beta','rotation',0); title(TT,'interpreter','latex');
set(gcf, 'Position',  [0.1*Lx, 0.2*Ly, 0.8*Lx, 0.6*Ly]);

%% Plot Stream Function
figure(2); TT = ['$\psi(\eta,\beta)$ Contour Plot: $\sqrt{R_e} =$ ',num2str(Re)];
contourf(eta,beta,Psi,15,'LineColor','none' ); 
colorbar('location','southoutside'); caxis([0,0.8]); colormap('jet');
xlabel('\eta'); ylabel('-\beta','rotation',0); title(TT,'interpreter','latex');
set(gcf, 'Position',  [0.1*Lx, 0.2*Ly, 0.8*Lx, 0.6*Ly]);

%% Plot Vorticity
figure(3); T = ['$\omega(\eta,\beta)$ Contour Plot: $\sqrt{R_e} =$ ',num2str(Re)];
contourf(eta,beta,Omega,15,'LineColor','none'); 
colorbar('location','southoutside'); caxis([-0.15,0.1]); colormap('jet');
xlabel('\eta'); ylabel('-\beta','rotation',0); title(T,'interpreter','latex');
set(gcf, 'Position',  [0.1*Lx, 0.2*Ly, 0.8*Lx, 0.6*Ly]);

%% Plot V velocity component 
figure(4); T = ['$V(\eta,\beta)$ at $\sqrt{R_e} =$ ',num2str(Re)];
fn = pcolor(eta,beta,V); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
colorbar('location','southoutside'); colormap('jet');
xlabel('\eta'); ylabel('-\beta','rotation',0); title(T,'interpreter','latex');
set(gcf, 'Position',  [0.1*Lx, 0.2*Ly, 0.8*Lx, 0.6*Ly])

%% Plot W velocity component 
figure(5); T = ['$W(\eta,\beta)$ at $\sqrt{R_e} =$ ',num2str(Re)];
Wp = W; Wp(Wp<0)=NaN; % remove negative/reverse flow
fn = pcolor(eta,beta,Wp); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
colorbar('location','southoutside'); caxis([0,0.15]); colormap('jet'); 
xlabel('\eta','fontsize',12); ylabel('-\beta','fontsize',12,'rotation',0); title(T,'interpreter','latex','fontsize',12);
set(gcf, 'Position',  [200, 200, 1200, 400])
%figfile = ['Stew_W_Re=',num2str(Re),'.png']; saveas(figure(5),figfile); pause(1)

%% Plot U velocity component 
figure(6); T = ['$U(\eta,\beta)$ at $\sqrt{R_e} =$ ',num2str(Re)];
fn = pcolor(eta,beta,-U); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
colorbar('location','southoutside'); caxis([0,0.15]); colormap('jet'); 
xlabel('\eta'); ylabel('-\beta','rotation',0); title(T,'interpreter','latex');
set(gcf, 'Position',  [0.1*Lx, 0.2*Ly, 0.8*Lx, 0.6*Ly])

%% Plot Pressure 
figure(7); T = ['$\bar{P}(\eta,\beta)$ at $\sqrt{R_e} =$ ',num2str(Re)];
fn = pcolor(eta,beta,P); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
colorbar('location','southoutside'); caxis([-0.015,0.01]); colormap('jet'); 
xlabel('\eta'); ylabel('-\beta','rotation',0); title(T,'interpreter','latex');
set(gcf, 'Position',  [0.1*Lx, 0.2*Ly, 0.8*Lx, 0.6*Ly])
