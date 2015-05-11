figure(3)
clear mm
%fix0(:,:)=Fm(Nx/2,:,:,1);
fl=size(Fm);
for ii=1:fl(4)
 ii  
% clf
%     fix(:,:)=Fm(Nx/2,:,:,ii);
%     contour(fix,[0 0],'k')
%     hold on
%     contour(fix0,[0 0],'r')
%         axis equal
%     hold off
%  
%     mm(ii)=getframe(gcf);
%     
  clf
fi(:,:,:)=Fm(:,:,:,ii);
fi = smooth3(fi,'box',5);
isosurface(fi,0)

view(28,36),
%camlight(10,10)
axis equal
axis([1 Nx 1 Ny 1 Nz])
mm(ii)=getframe(gcf);
 end
movie2avi(mm,'2rename','compression','none')