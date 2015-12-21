% Program to calculate turing BVAM 3 dimensions

clear all

 NF=200;
 eta=sqrt(2);
 Du=.516/eta;
 Nx=30;
 Ny=30;
 Nz=30;
 R=11;

h=-1.;
a=1/.899;
b=-.91/.899;
c=0.02;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% form of fi and initial values
%load dominio
fiini=fi;
%%%%%%%%%%% parameters for iteraion loop %%%%%%%%%%%%%%%%%%%%%%%%
step=400;
iter=1;
dt=.02;
cont=iter;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initial conditions for u and beta %%%%%
        u=.2*(rand(Nx,Ny,Nz)-.5);
        v=.2*(rand(Nx,Ny,Nz)-.5);        
% 
         u(fi<=-.9)=0;
         v(fi<=-.9)=0;
%%    store the initial domain

fix0(:,:)=fi(Nx/2,:,:);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

%%
load turing2
for iter=cont:NF              %time loop
    for iiter=1:step
 
 %%   Laplacoano con |nabla^2=oen la frontera del cubo 
        H=u;
        lap3Dt
        lapu=lapH;
        

        H=v;
        lap3Dt
        lapv=lapH; 
        
        lapu(fi<=-.9)=0;
        lapv(fi<=-.9)=0;
        u(fi<=-.9)=0;
        v(fi<=-.9)=0;

%       u=u+dt*(Du*lapu+(u+a*v-c*u.*v-u.*v.^2));

        u=u+dt*(Du*(lapu+lapFu)+S+(u+a*v-c*u.*v-u.*v.^2));
        v=v+dt*(lapv/eta+(b*v+h*u+c*u.*v+u.*v.^2));        
%        u(:,:,1)=u(:,:,2);
%       v(:,:,1)=v(:,:,2);
    end
    ux(:,:)=u(:,Ny/2,:);
    fix(:,:)=fi(Nx/2,:,:);

    %%   
hh=max(max(max(isnan(fi(:,:,:)))));
    if hh==1;
       'nans'
        break
    end
      
    iter
    
    %%
                
figure(1)
clf
%fim = smooth3(fi,'box',5);
%isosurface(S)
%camlight(10,40,'infinity')

% axis equal
% axis([1 Nx 1 Ny 1 Nz])
% light
% material metal
% colorbar
hold on
surf(ux,'FaceAlpha',0.5),shading interp,  view(-36,18)
%mesh(fix,'FaceColor','none')
hold off
 %%
      figure(2)
     clf
[x,y,z] = meshgrid(1:1:Ny,1:1:Nx,1:1:Nz);
xslice = [Nx/2-R:R:Nx/2+R,Nx/2-R:R:Nx/2+R];yslice = [Ny/2-R:R:Ny/2+R,Ny/2-R:R:Ny/2+R]; zslice = [0:3:R,0:3:R];
p3=slice(x,y,z,u,xslice,yslice,zslice);
set(p3,'FaceColor','flat','EdgeColor','none','FaceAlpha',0.1);
rs=max(abs(max(max(max(u)))),abs(min(min(min(u)))));

axis equal, view(74,18), 
set(gca,'CLim',[-rs,rs])
colorbar;
%%
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
    light
material metal
    colorbar
view(-15,40)


%para hacer peliculas 3D usar cine3D y cine2D

    Fm(:,:,:,iter)=fi(:,:,:);
    U(:,:,:,iter)=u(:,:,:);
    Sm(:,:,:,iter)=S(:,:,:);
    %%
%     figure(4)
%     fix(:,:)=fi(:,Ny/2,:);
%     contour(fix,[0 0],'k')
%     hold on
%     contour(fix0,[0 0],'r')
%     axis equal
%     %getframe(gcf);
%     hold off
    
end

