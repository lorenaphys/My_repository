load('julio15h.mat');

f = size(Um);
Nx = f(1);
Ny = f(2);
Nz = f(3);
h = 101;
fi0 = Fm(:,:,:,1);
Fi = 0*Um;
for i = 1:cont+1
	Fi(:,:,:,i) = fi0;
end
for j = cont+2:cont+2+NF 
	Fi(:,:,:,j) = Fm(:,:,:,j-cont-1);
end
for k = cont+2+NF:cont+1+NF+cont3
	Fi(:,:,:,k) = Fm(:,:,:,NF+1);
end
%Fi(:,:,:,cont+2:cont+2+NF) = Fm(:,:,:,2:NF+1);
%Fi(:,:,:,cont+2+NF:cont+1+NF+cont3) = Fm(:,:,:,NF+1);
%M = struct('cdata',[],'colormap',[]);
M1 = struct('cdata',[],'colormap',[]);

for i = 1:f(4)
     clf
     cdata = smooth3((Um(:,:,:,i)-min(min(min(Um(:,:,:,i)))))./...
             (max(max(max(Um(:,:,:,i))))-min(min(min(Um(:,:,:,i))))),'box',5);
     fi = smooth3(Fi(:,:,:,i),'box',5);
     p4=patch(isosurface(fi,0));
     isonormals(fi,p4);
     isocolors(cdata,p4);
     set(p4,'FaceColor','interp','EdgeColor','none'),
     camlight, lighting phong
     axis equal, view(-16,24), axis off 
     axis([1 40 1 40 1 25]),
     colormap jet
     M1(i) = getframe;
     disp(i)
end

% for k = 101:f(4)
%     R = 11;
%     u = Um(:,:,:,k);
% %     [x,y,z] = meshgrid(1:1:Ny,1:1:Nx,1:1:Nz);
% %     xslice = [Nx/2-R:R:Nx/2+R,Nx/2-R:R:Nx/2+R];yslice = [Ny/2-R:R:Ny/2+R,Ny/2-R:R:Ny/2+R]; zslice = [0:3:R,0:3:R];
% %     p3=slice(x,y,z,u,xslice,yslice,zslice);
% %     set(p3,'FaceColor','flat','EdgeColor','none','FaceAlpha',0.1);
% %     rs=max(abs(max(max(max(u)))),abs(min(min(min(u)))));
% %     axis equal, view(74,18), 
% %     set(gca,'CLim',[-rs,rs])
% %     colormap hsv;
% 
%     cdata = smooth3((Um(:,:,:,k)-min(min(min(Um(:,:,:,k)))))./...
%             (max(max(max(Um(:,:,:,k))))-min(min(min(Um(:,:,:,k))))),'box',5);
%     [x,y,z] = meshgrid(1:1:Nx,1:1:Ny,1:1:Nz);
%     xslice = [0,Nx/2,Nx];yslice = [0,Ny/2,Ny]; zslice = [0,Nz/2,Nz];
%     p3=slice(x,y,z,Um(:,:,:,k),xslice,yslice,zslice);
%     set(p3,'FaceColor','interp','EdgeColor','none','FaceAlpha',0.5),
%     axis equal, view(-70,20)
%     colormap jet,
%     M(k) = getframe;
%     disp(k)
% end
