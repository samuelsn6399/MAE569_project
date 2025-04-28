function [] = plot_orbit(r,v,e,a,Omega,omega,inc,isdeg,unit)
%plot_orbit plots the orbit of a spacecraft given that the orbit is
%elliptical.
if nargin == 8; unit = 'DU'; end

xr = [0,r(1)];
yr = [0,r(2)];
zr = [0,r(3)];
xv = [r(1),r(1)+v(1)];
yv = [r(2),r(2)+v(2)];
zv = [r(3),r(3)+v(3)];
transform_ellipse = orbit_ellipse_coords(a,e,Omega,omega,inc,isdeg);

limit = max(abs([transform_ellipse(1,:);transform_ellipse(2,:)]),[],'all');

figure
hold on
surf([-limit,-limit;limit,limit],[-limit,limit;-limit,limit],zeros([2,2]),'FaceAlpha',0.3,FaceColor='Blue')
plot3(xr,yr,zr,Color='red',LineWidth=2)
plot3(r(1),r(2),r(3),'^r',MarkerSize=5,LineWidth=2)
plot3(xv,yv,zv,Color='red',LineWidth=2)
plot3(r(1)+v(1),r(2)+v(2),r(3)+v(3),'^r',MarkerSize=5,LineWidth=2)
plot3(transform_ellipse(1,:),transform_ellipse(2,:),transform_ellipse(3,:),'r-',LineWidth=2)
plot3([0,5],[0,0],[0,0],'k-',LineWidth=1.5)
plot3([0,0],[0,5],[0,0],'k-',LineWidth=1.5)
plot3([0,0],[0,0],[0,5],'k-',LineWidth=1.5)
view([53.848 44.751])
[xe,ye,ze] = sphere;
if unit == 'km'
    xe = xe*6378.136;
    ye = ye*6378.136;
    ze = ze*6378.136;
elseif unit == 'mi'
    xe = xe*3963.19;
    ye = ye*3963.19;
    ze = ze*3963.19;
end
surf(xe,ye,ze,FaceColor='Green')
axis equal
hold off
view([1 1 1])
end

