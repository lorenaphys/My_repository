
load('feb3c.mat')

fl=size(Fm);
Nx = fl(1);
Ny = fl(2);
Nz = fl(3);
f = 66;
fi(:,:,:)=Fm(:,:,:,f);
u(:,:,:)=Um(:,:,:,f);   
clf
       
figure(1)   
     cdata = smooth3((u-min(min(min(u))))./(max(max(max(u)))-min(min(min(u)))),'box',5);
     fi = smooth3(fi,'box',5);
     p4=patch(isosurface(fi,0));
     isonormals(fi,p4);
     isocolors(cdata,p4);
     set(p4,'FaceColor','interp','EdgeColor','none'),
     camlight, lighting phong
     axis equal, view(-14,40), axis off
     %axis([4.5 35 4.5 35 1 19]),
     axis([1 Nx 1 Ny 1 Nz])
     colormap jet
     %colorbar
     

figure(2)
   fix0(:,:) = Fm(:,Ny/2,:,1);
   fix(:,:)=fi(:,Ny/2,:);
   contour(fix,[0 0],'k')
   hold on
   contour(fix0,[0 0],'r')
   axis equal
   axis([1 20 3 36])
    %getframe(gcf);
%    hold off
    
% clf
% fi = smooth3(fi,'box',5);
% isosurface(fi,0)
% view(64,26),
% axis equal
% axis([1 Nx 1 Ny 1 Nz])
% mm(:,:,ii)=getframe(gcf);

% clf
% cdata = smooth3((u-min(min(min(u))))./(max(max(max(u)))-min(min(min(u)))),'box',5);
% [x,y,z] = meshgrid(1:1:Nx,1:1:Ny,1:1:Nz);
% xslice = [Nx/2,Ny/2];yslice = Ny/2; zslice = [0,10];
% p3=slice(x,y,z,u,xslice,yslice,zslice);
% set(p3,'FaceColor','interp','EdgeColor','none','FaceAlpha',0.5),
% axis equal, view(-70,20)
%  mm(:,:,ii)=getframe(gcf); 

%exit   

% figure(3)
% t = linspace(1,fl(4),fl(4));
% plot(Sm)
% xlabel('t')
% ylabel('\sigma','FontSize',14,'FontWeight','bold')
% axis([0 fl(4) 0 25])

u1 =  Um(:,:,Nz/2,1);
figure(4)
%[X,Y] = meshgrid(1:Nx,1:Ny);
surf(u1)

uf = u(:,:,Nz/2);
figure(5)
surf(uf)

