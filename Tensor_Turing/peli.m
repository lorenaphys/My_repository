load('junio10d.mat');

f = size(Fm);
Nx = f(1);
Ny = f(2);
Nz = f(3);
h = 51;
M1 = struct('cdata',[],'colormap',[]);

for i = 1:h
     clf
     cdata = smooth3((Um(:,:,:,i)-min(min(min(Um(:,:,:,i)))))./...
             (max(max(max(Um(:,:,:,i))))-min(min(min(Um(:,:,:,i))))),'box',5);
     fi = smooth3(Fm(:,:,:,i),'box',5);
     p4=patch(isosurface(fi,0));
     isonormals(fi,p4);
     isocolors(cdata,p4);
     set(p4,'FaceColor','interp','EdgeColor','none'),
     camlight, lighting phong
     axis equal, view(-16,24), axis off 
     axis([1 40 1 40 1 70]),
     colormap jet
     M1(i) = getframe;
end