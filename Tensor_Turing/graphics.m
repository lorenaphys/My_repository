
load('junio27d.mat')


fl=size(Fm);
Nx = fl(1);
Ny = fl(2);
Nz = fl(3);
f = 7;
M = struct('cdata',[],'colormap',[]);
fi0(:,:,:) = Fm(:,:,:,1);
fi(:,:,:)=Fm(:,:,:,f);
u(:,:,:)=Um(:,:,:,f);
u0(:,:,:) = Um(:,:,:,1);
%u0(fi0<=-0.99) = 0;
clf

%%
figure(1)   
     cdata = smooth3((u0-min(min(min(u0))))./(max(max(max(u0)))-min(min(min(u0)))),'box',5);
     fi0 = smooth3(fi0,'box',5);
     p4=patch(isosurface(fi0,0));
     isonormals(fi0,p4);
     isocolors(cdata,p4);
     set(p4,'FaceColor','interp','EdgeColor','none'),
     camlight, lighting phong
     axis equal, view(-16,24), axis off %de (18,10) lo cambie a (-16,24)%segundo cambio para view, antes era (-16,24)%cambié view de (-14,40) a (-16,24)
     axis([1 40 1 40 1 70]),
     %axis([1 Nx 1 Ny 1 Nz])
     colormap jet
     %colorbar

%%
figure(2)   
     cdata = smooth3((u-min(min(min(u))))./(max(max(max(u)))-min(min(min(u)))),'box',5);
     fi = smooth3(fi,'box',5);
     p4=patch(isosurface(fi,0));
     isonormals(fi,p4);
     isocolors(cdata,p4);
     set(p4,'FaceColor','interp','EdgeColor','none'),
     camlight, lighting phong
     axis equal, view(-16,24), axis off %de (18,10) lo cambie a (-16,24)%segundo cambio para view, antes era (-16,24)%cambié view de (-14,40) a (-16,24)
     axis([1 40 1 40 1 70]),
     %axis([1 Nx 1 Ny 1 Nz])
     colormap jet
     %colorbar
     
%%
figure(3)
   fix0(:,:) = Fm(:,Ny/2,:,1);
   fix(:,:)=fi(:,Ny/2,:);
   contour(fix,[0 0],'k')
   hold on
   contour(fix0,[0 0],'r')
   axis equal
   axis([1 60 1 40])
   legend('\phi','\phi_0'),
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
% axis([0 fl(4) 0 100])

%%
u1 =  Um(:,:,Nz/2,1);
figure(4)
%[X,Y] = meshgrid(1:Nx,1:Ny);
surf(u1)

%%
uf = u(:,:,Nz/2);
figure(5)
surf(uf)

%%
% figure(6)
% for k = 1:fl(4)
%     
%     cdata = smooth3((Um(:,:,:,k)-min(min(min(Um(:,:,:,k)))))./...
%            (max(max(max(Um(:,:,:,k))))-min(min(min(Um(:,:,:,k))))),'box',5);
%     [x,y,z] = meshgrid(1:1:Nx,1:1:Ny,1:1:Nz);
%     xslice = [0,10,20,30];yslice = [0,10,20,30]; zslice = [20,30,40];
%     p3=slice(x,y,z,Um(:,:,:,k),xslice,yslice,zslice);
%     set(p3,'FaceColor','interp','EdgeColor','none','FaceAlpha',0.5),
%     axis equal, view(-70,20)
%     colormap jet,
%     M(k) = getframe;
%     disp(k)
% end
