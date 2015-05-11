figure(3)
clear mm
%fix0(:,:)=Fm(Nx/2,:,:,1);
fl=size(Fm);
for ii=1:fl(4)
ii
fiaux(:,:,:)=Fm(:,:,:,ii);
betaux(:,:,:)=Um(:,:,:,ii);
%%
figure(1)
clf
%fiaux = smooth3(fiaux,'box',5);
isosurface(fiaux,0)
view(64,26),
axis equal
axis([1 Nx 1 Ny 1 Nz])


%%
    figure(2)
%     clf
%     cdata = smooth3((betaux-min(min(min(betaux))))./(max(max(max(betaux)))-min(min(min(betaux)))),'box',5);
%      betaux = smooth3(betaux,'box',5);
% 
% cdata = smooth3((betaux-min(min(min(betaux))))./(max(max(max(betaux)))-min(min(min(betaux)))),'box',5);
% [x,y,z] = meshgrid(1:1:Nx,1:1:Ny,1:1:Nz);
% xslice = [Nx/2,Nx/2];yslice = [Ny/2-10,Ny/2+10]; zslice = [R,2*R];
p3=slice(x,y,z,betaux,xslice,yslice,zslice);
set(p3,'FaceColor','interp','EdgeColor','none','FaceAlpha',0.5);
axis equal, view(90,90) %view(68,16), 
colorbar
pause(.1)
% 
%           M(:,:,iter)=getframe;


    %%
   figure(3)
   fix(:,:)=fiaux(:,Ny/2,:);
   contour(fix,[0 0],'k')
   hold on
   contour(fix0,[0 0],'r')
   axis equal
    getframe(gcf);
   hold off
    
    
    
    

% clf
%     fix(:,:)=Fm(Nx/2,:,:,ii);
%     contour(fix,[0 0],'k')
%     hold on
%     contour(fix0,[0 0],'r')
%         axis equal
%     hold off
%  
%     mm(ii)=getframe(gcf);
    
%   clf
% fi(:,:,:)=Fm(:,:,:,ii);
% fi = smooth3(fi,'box',5);
% isosurface(fi,0)
% 
% view(28,36),
% %camlight(10,10)
% axis equal
% axis([1 Nx 1 Ny 1 Nz])
% mm(ii)=getframe(gcf);
 end
%movie2avi(mm,'2rename','compression','none')