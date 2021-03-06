Nx=40;
Ny=40;
Nz=70;
R=15;

% load initial
% u0=uout+ues; v0=vout+ves; w0=wout+wes; x0=xout+xes;
fi=ones(Nx,Ny,Nz);
 r = zeros(Nx,Ny,Nz);
% 
% %ancho=2;
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            r(i,j,k)=sqrt((i-Nx/2)^2+(j-Ny/2)^2+(k)^2);
            if r(i,j,k)>=R
            fi(i,j,k)=-1;
            end
        end
   
   end
end
load('oct25e.mat')
u = Um(:,:,:,301);

    %R = 9;
    %N = 0;
    %[X,Y,Z]=meshgrid(1:Nx,1:Ny,1:Nz);
        %teta=atan2((Y-Ny/2),(X-Nx/2));
        %rad=sqrt((X-Nx/2+.5).^2+(Y-Ny/2+.5).^2);
        %u=1.5*exp(-((X-Nx/3-.5).^2+(Y-Ny/3-.5).^2+(Z-13.5).^2)/50);
        %u=2.5*rad.*(cos(teta*N)+sin(teta*N)).*(Z/Nz)/max(max(max(rad)))+(exp(-((-X+Nx/2-.5).^2+(-Y+Ny/2-.5).^2+(-Z+R+14).^2)/80));
        %u=-2.5*rad.*(cos(teta*N)+sin(teta*N)).*(Z/Nz)/max(max(max(rad)))+(exp(-((X-Nx/2).^2+(Y-Ny/2).^2+(Z-R+14).^2)/50));
        %u=2.5*rand(Nx,Ny,Nz);
        %u = 0.2*(rand(Nx,Ny,Nz)-.5);
        %u=.1*u+.2*(rand(Nx,Ny,Nz)-.5);        

% fix0(:,:)=fi(Nx/2,:,:);
% 
% figure(1)
% clf
%     ux(:,:)=u(:,Ny/2,:);
%     %vx(:,:)=v(:,Ny/2,:);
%     fix(:,:)=fi(Nx/2,:,:);
% clf
% hold on
% mesh(ux,'FaceAlpha',0.5),shading interp,  view(-36,18)
% %mesh(fix,'FaceColor','none')
% hold off
% 
% % %%
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
    figure(1)
    
%     u = smooth3(u,'box',3); 
    clf
       
    cdata = smooth3((u-min(min(min(u))))./(max(max(max(u)))-min(min(min(u)))),'box',5);
       %cdata = smooth3(u,'box',3);
    fim = smooth3(fi,'box',3);
    p4=patch(isosurface(fim,0));
    isonormals(fim,p4);
    isocolors(cdata,p4);
    set(p4,'FaceColor','interp','EdgeColor','none'),
    camlight, lighting phong
    axis equal, axis off, 
    axis([1 40 1 40 1 20]),
    colormap jet
    %light
%material metal
    %colorbar
%view(-15,40)
view(-16,24)


%     figure(4)
%     fix(:,:)=fi(:,Ny/2,:);
%     contour(fix,[0 0],'k')
%     hold on
%     contour(fix0,[0 0],'r')
%     axis equal
%     %getframe(gcf);
%     hold off
