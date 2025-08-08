clear;close all;clc
n=2000;
[r,theta,phi]=generatemesh(1,n);

ax=figure('Color','White');
pos=[1 1 15 15];
set(gcf,'unit','centimeters','position',pos)
ax1=axes('position',[0,0.5,0.5,0.5])  ;
sphere;
% cmap=hot(256);
% cmap=cmap(32:96,:);
colormap("gray")
axis equal
axis off
hold on
x=r.*sin(theta).*cos(phi);
y=r.*sin(theta).*sin(phi);
z=r.*cos(theta);
text(0-0.1,0,-1.4,'(a)',"FontSize",13,'FontWeight','bold')
% shading flat
plot3(x,y,z,'.','Color','r')%,'Color','b'
ax1=axes('position',[0.5,0.6,0.35,0.3])  ;
r1=r(1:n/2).*sin(theta(1:n/2));
%r1=r1(r1)
polarplot(phi(1:n/2),r1,'.','Color','r');
set(gca,'ThetaTickLabel',[]);
set(gca,'RTickLabel',[]);
set(gca,'RTick',[]);
set(gca,'GridColor','k');
set(gca,'GridAlpha',0.45);
axis on
hold on
text(-pi/2-0.1,1.315,'(b)',"FontSize",13,'FontWeight','bold')
%polarplot(linspace(0,2*pi,n/2),ones(n/2));
%export_fig ESDmesh.png -r600