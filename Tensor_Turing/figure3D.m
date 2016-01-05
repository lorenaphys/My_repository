load('prueba.mat');

f = size(Fm);
Nx = f(1);
Ny = f(2);
Nz = f(3);
fl = f(4);
R = 15;

fi(:,:,:) = Fm(:,:,:,fl);
u(:,:,:) = Um(:,:,:,fl);

fix0(:,:)=fi(Nx/2,:,:);

figure(1)
clf
    ux(:,:)=u(:,Ny/2,:);
    %vx(:,:)=v(:,Ny/2,:);
    fix(:,:)=fi(Nx/2,:,:);
clf
hold on
mesh(ux,'FaceAlpha',0.5),shading interp,  view(-36,18)
%mesh(fix,'FaceColor','none')
hold off

% %%
%       figure(2)
%      clf
% [x,y,z] = meshgrid(1:1:Ny,1:1:Nx,1:1:Nz);
% xslice = [Nx/2-R:R:Nx/2+R,Nx/2-R:R:Nx/2+R];yslice = [Ny/2-R:R:Ny/2+R,Ny/2-R:R:Ny/2+R]; zslice = [0:3:R,0:3:R];
% p3=slice(x,y,z,u,xslice,yslice,zslice);
% set(p3,'FaceColor','flat','EdgeColor','none','FaceAlpha',0.1);
% rs=max(abs(max(max(max(u)))),abs(min(min(min(u)))));
% 
% axis equal, view(74,18), 
% set(gca,'CLim',[-rs,rs])
% colorbar;
% %%
    figure(3)
    
%     u = smooth3(u,'box',3); 
    clf
       
    %cdata = smooth3((u-min(min(min(u))))./(max(max(max(u)))-min(min(min(u)))),'box',5);
       cdata = smooth3(u,'box',3);
    fim = smooth3(fi,'box',3);
    p4=patch(isosurface(fim,0));
    isonormals(fim,p4);
    isocolors(cdata,p4);
    set(p4,'FaceColor','interp','EdgeColor','none'),
    camlight, lighting phong
    axis equal, axis off, 
    axis([1 Nx 1 Ny 1 Nz]),
    colormap jet
    light
material metal
    colorbar
view(-15,40)


    figure(4)
    fix(:,:)=fi(:,Ny/2,:);
    contour(fix,[0 0],'k')
    hold on
    contour(fix0,[0 0],'r')
    axis equal
    %getframe(gcf);
    hold off
