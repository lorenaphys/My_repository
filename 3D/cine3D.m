clear mm
fl=size(Fm);
for ii=1:fl(4)
 ii  
%      hold on 
%      Xn=frame2im(M(ii));
%      image(Xn)

 fi(:,:,:)=Fm(:,:,:,ii);
  bet(:,:)=U(:,:,ii);   
    clf
       
    cdata = smooth3((bet-min(min(min(bet))))./(max(max(max(bet)))-min(min(min(bet)))),'box',5);
    fi = smooth3(fi,'box',5);
    p4=patch(isosurface(fi,0));
    isonormals(fi,p4);
    isocolors(cdata,p4);
    set(p4,'FaceColor','interp','EdgeColor','none'),
    camlight, lighting phong
    axis equal, view(-14,40), axis off
    axis([1 Nx 1 Ny 1 Nz/2]),
    colorbar
    mm(:,:,ii)=getframe(gcf);

% 
% clf
% fi = smooth3(fi,'box',5);
% isosurface(fi,0)
% view(64,26),
% axis equal
% axis([1 Nx 1 Ny 1 Nz])
% mm(:,:,ii)=getframe(gcf);

% clf
% cdata = smooth3((bet-min(min(min(bet))))./(max(max(max(bet)))-min(min(min(bet)))),'box',5);
% [x,y,z] = meshgrid(1:1:Nx,1:1:Ny,1:1:Nz);
% xslice = [Nx/2,Ny/2];yslice = Ny/2; zslice = [0,10];
% p3=slice(x,y,z,bet,xslice,yslice,zslice);
% set(p3,'FaceColor','interp','EdgeColor','none','FaceAlpha',0.5),
% axis equal, view(-70,20)
%  mm(:,:,ii)=getframe(gcf);    
end

movie2avi(mm,'1rename','compression','none')