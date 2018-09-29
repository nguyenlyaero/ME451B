clear, clc
close all
global sigma s r

%% Input parameters
sigma=10; s=8/3; % classical case
r=0.1;          % Try: r = 0.1, r = 1.24, r = 10, r = 22, r = 26

% Observations with increasing r (for fixed [sigma = 10, s = 8/3]): 
% 1. Globally stable for 0<r<1
% 2. Pitchfork bifurcation at r=1
% 3. Node-spiral transition at r=1.246 (oscillatory behavior r>1.246)
% 4. Homoclinic Orbits for r>13.93 (solution loops between the two stable fixed points)
% 5. Hopf Bifurcation at r=24.06 (chaotic solution for r>24.06)
% 6. What happens after r>24.06? (Why does the solution stay bounded?)
purt_mag = 0.01; % purturbation magnitude
r_Hopf = sigma*(sigma + s + 3)/(sigma - s - 1) % When does the Hopf Bifurcation occur? 

%% Initialize (purturbation)
Purt = purt_mag*rand(3,1); 
xpurt = Purt(1); ypurt = Purt(2); zpurt = Purt(3); % Purturbation magnitude
fixed_point_init = 1; % perturb about which fixed point? 1: Null fixed point
tend = 100;

%% Get the fixed points and plot them
[xfixed,yfixed,zfixed] = FixedPointsLorenz_3eq(s,r);
figure(1), hold on
for idx = 1:length(xfixed)
    figure(1), plot3(xfixed(idx),yfixed(idx),zfixed(idx),'ok','MarkerFaceColor','k')
end

%% Solve the dynamical system
X0 = xfixed(fixed_point_init) + xpurt; 
Y0 = yfixed(fixed_point_init) + ypurt; 
Z0 = zfixed(fixed_point_init) + zpurt;
[t,XYZ] = ode45(@RHSLorenz_3eq,[0 tend],[X0;Y0;Z0]);
X = XYZ(:,1); Y = XYZ(:,2); Z = XYZ(:,3);

%% Plot the solution
figure(1), plot3(X,Y,Z,'linewidth',0.5)
view([1 -1 1.2]), xlabel('x'), ylabel('y'), zlabel('z')
plot3(X(1),Y(1),Z(1),'*g')
plot3(X(end),Y(end),Z(end),'*r')
set(gca,'fontsize',16), box on, daspect([1 1 1])
xlim([min(xfixed) - r, max(xfixed) + r])
ylim([min(yfixed) - r, max(yfixed) + r])
zlim([min(zfixed) - r, max(zfixed) + r])
%view([0 0 1])
if (r > 1) 
    figure(1), legend('Fixed Pt. 1', 'Fixed Pt. 2', 'Fixed Pt. 3','Solution Trajectory','Initial Location','Final Location')
else
    figure(1), legend('Fixed Pt. 1','Solution Trajectory','Initial Location','Final Location')
end
    
%% Check the stability of the fixed points
for idx = 1:length(xfixed)
    e = JacobianLorenz_3eq(sigma,r,s,xfixed(idx),yfixed(idx),zfixed(idx));
    display(strcat('Eigenvalue at Fixed point #',num2str(idx)));
    display(e)
end

%% Plot the time history of the coordinates
figure(2), subplot(3,1,1), plot(t,X), xlabel('t'), ylabel('X')
set(gca,'fontsize',16), box on
figure(2), subplot(3,1,2), plot(t,Y), xlabel('t'), ylabel('Y')
set(gca,'fontsize',16), box on
figure(2), subplot(3,1,3), plot(t,Z), xlabel('t'), ylabel('Z')
set(gca,'fontsize',16), box on
