load('oct25i.mat');

f = size(Um);
Nx = f(1);
Ny = f(2);
Nz = f(3);
h = 101;
%fi0 = Fm(:,:,:,1);
%Fi = 0*Um;
%for i = cont1:cont2
	%Fi(:,:,:,i) = fi0;
%end
%for j = cont+2:cont+2+NF 
	%Fi(:,:,:,j) = Fm(:,:,:,j-cont-1);
%end
%for k = cont2+cont4+2:cont2+cont4+2+cont5
	%Fi(:,:,:,k) = Fm(:,:,:,cont2+cont4+1);
%end
M = struct('cdata',[],'colormap',[]);
%M1 = struct('cdata',[],'colormap',[]);

%for i = 1:f(4)
     %clf
     %cdata = smooth3((Um(:,:,:,i)-min(min(min(Um(:,:,:,i)))))./...
             %(max(max(max(Um(:,:,:,i))))-min(min(min(Um(:,:,:,i))))),'box',5);
     %fi = smooth3(Fm(:,:,:,i),'box',5);
     %p4=patch(isosurface(fi,0));
     %isonormals(fi,p4);
     %isocolors(cdata,p4);
     %set(p4,'FaceColor','interp','EdgeColor','none'),
     %camlight, lighting phong
     %axis equal, view(-16,24), axis off 
     %axis([1 40 1 40 1 30]),
     %colormap jet
     %M1(i) = getframe;
     %disp(i)
%end

%for k = 1:f(4)
   %R = 11;
   %u = Um(:,:,:,k); 
    %cdata = smooth3((Um(:,:,:,k)-min(min(min(Um(:,:,:,k)))))./...
            %(max(max(max(Um(:,:,:,k))))-min(min(min(Um(:,:,:,k))))),'box',5);
    %[x,y,z] = meshgrid(1:1:Nx,1:1:Ny,1:1:Nz);
    %xslice = [0,Nx/2];yslice = [0,Ny/2]; zslice = [0,10];
    %p3=slice(x,y,z,Um(:,:,:,k),xslice,yslice,zslice);
    %set(p3,'FaceColor','interp','EdgeColor','none','FaceAlpha',0.5),
    %axis equal, view(-70,20)
    %colormap jet,
    %M(k) = getframe;
    %disp(k)
%end

for k = 301:301
    clf
    u = smooth3(Um(:,:,:,k),'box',5);
	isosurface(u,(u-min(min(min(u))))./...
            (max(max(max(u)))-min(min(min(u)))));
	axis off, axis equal, view(-16,24)
	axis([1 50 1 50 1 70])
	colormap winter
	%camlight
	M(k) = getframe;
	disp(k)
end
