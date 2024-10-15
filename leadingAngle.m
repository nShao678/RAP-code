clear;
clc;
close all
figure;
hold on;
axis off;
% grid on;



theta = linspace(0, 2*pi, 100);
plot(cos(theta), sin(theta), 'k--', 'LineWidth', 1,'HandleVisibility','off');
theta0 = pi/3;
dtheta = pi/4;
theta = theta0+linspace(-dtheta,dtheta, 100);
plot(cos(theta), sin(theta), 'k-', 'LineWidth', 2,'DisplayName','$\mathsf{Rq}(x)\leq \rho$');
plot(cos(theta0),sin(theta0),'ko', 'LineWidth', 2,'DisplayName','$u_{1}$');
plot(cos(theta0+dtheta/2),sin(theta0+dtheta/2),'ro', 'LineWidth', 2,'DisplayName','$x$');
t = linspace(-2, 2, 100);
plot(cos(theta0+dtheta/2)-t*sin(theta0+dtheta/2),sin(theta0+dtheta/2)+t*cos(theta0+dtheta/2),'r-', 'LineWidth', 2,'DisplayName','$\{v| v^{\mathsf{T}}Ax=0\}$');
plot(cos(theta0+dtheta/2)-t*sin(theta0+dtheta/6),sin(theta0+dtheta/2)+t*cos(theta0+dtheta/6),'b-', 'LineWidth', 2,'DisplayName','$\{v| v^{\mathsf{T}}Bx=0\}$');
t = linspace(-0.3, 2, 100);
plot(cos(theta0+dtheta/2)+t*cos(theta0+dtheta/2),sin(theta0+dtheta/2)+t*sin(theta0+dtheta/2),'r--', 'LineWidth', 1,'HandleVisibility','off');

arc_theta = linspace(pi/2+(theta0+dtheta/6),theta0+dtheta/2, 50);
arc_radius = 0.3;
arc_x = arc_radius * cos(arc_theta);
arc_y = arc_radius * sin(arc_theta);
plot(arc_x+cos(theta0+dtheta/2), arc_y+sin(theta0+dtheta/2), 'g', 'LineWidth', 3,'DisplayName','Leading angle $\vartheta$');

hold off
axis([-0.5,1,0,1.4])
set(gcf, 'Color', 'w');
legend('Interpreter','latex','FontSize',14,'Location','southwest','Box','off')

export_fig ('leadingAngle.eps')
export_fig ('leadingAngle.pdf')