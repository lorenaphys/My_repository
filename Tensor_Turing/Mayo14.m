% Program to calculate phase fiels in 3 dimensions

clear all

et=.01;
dx=1;
NF=200;
sig=0*(1:NF);
ep1=2;
ep=ep1^2;
sigma=.1;
Nx=30;
Ny=30;
Nz=30;
R=11;
N=0;
sifiu=0.;
duu=.5; 
Du=50; 


%%%%%%%%%%%%%%%%%%%%% strength of the fields  %%%%%%%%%%%
Aout=2; 
Av=2; 
Afi=.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% form of fi and initial values

semiesf3D
%formas3D
fiini=fi;
%%%%%%%%%%% parameters for iteraion loop %%%%%%%%%%%%%%%%%%%%%%%%
step=400;
iter=1;
dt=1e-5;
cont=iter;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initial conditions for u and beta %%%%%

u1=1; 
u2=-.0; 

uout=0.0; 

bet=.5;

% uinit=0*fi;
% [uf vf]=find(fi>=-.9);
% for i=1:length(uf)
%  uinit(uf(i),vf(i))=u0+.1*(rand-.5);%+u(uf(i),vf(i));
% end
[X,Y,Z]=meshgrid(1:Nx,1:Ny,1:Nz);
        teta=atan2((Y-Ny/2),(X-Nx/2));
        rad=sqrt((X-Nx/2+.5).^2+(Y-Ny/2+.5).^2);
        u=3.5*+exp(-((X-Nx/2-.5).^2+(Y-Ny/2-.5).^2+(Z-R+2).^2)/20);
        u=1-1.5*(u);
%u=1.5*u0;
        u(fi<=-.9)=0;
        u0=sum(sum(sum(u)))/Nx/Ny/Nz;
 %       u(fi>=.9)=0;        
%u=smooth3(uinit,'gaussian',3);
%%    store the initial domain

%fi = smooth3(fi,'box',5);
fix0(:,:)=fi(Nx/2,:,:);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

%%
%load phase3d-2
for iter=cont:NF
    for iiter=1:step
       % u=0*u;
                
        H=fi;
        lap3D
        lapfi=lapH;

        r1=(u-u1).^2;
        r2=(u-u2).^2; 
        B=r1.*r2;
        rout=(u-uout).^2;
        
        mu=((fi-ep1.*(bet*u.^2)).*((fi).^2-1)-ep*lapfi);        
        
        H=mu;
        lap3D
        lapmu=lapH;
        
        F=Afi*((3*fi.^2-1-2*ep1*fi.*(bet*u.^2)).*mu-ep*lapmu);
        F=F+rout*Aout.*(fi);
        F=F+2*Av*fi.*(fi.^2-1).*B;        
        
        
%         H=F;
%         lap3D
%         lapF=lapH;

        H=u;
        lap3D
        lapu=lapH;
        
        Fs=lapfi;
        
        
        H=Fs;
        lap3D
        lapFs=lapH;
         %% seccion para el parametro de lagrange
%         H=lapF;
%         grad3D
%         derilapF=g0iH;
%         derjlapF=g0jH;
%         derklapF=g0kH;
%         
%         H=lapFs;
%         grad3D
%         derilapFs=g0iH;
%         derjlapFs=g0jH;
%         derklapFs=g0kH;
%         
%         H=fi;
%         grad3D
%         derifi=g0iH;
%         derjfi=g0jH;
%         derklfi=g0kH; 
%         
%         I=sum(sum(derifi.*derilapF+derjfi.*derjlapF+derkfi.*derklapF));
%         Is=sum(sum(derifi.*derilapFs.*rr+derjfi.*derjlapFs.*rr));
        
%         sigma=-I/Is;
%%
%%%   dinamica de la u  %%%%%%%%%%%%%

       Gu=Av*(fi.^2-1).^2.*((u-u1).*r2+(u-u2).*r1);
       Gu=Gu-Afi*ep1*bet*u.*mu.*((fi).^2-1); 
       Gu=Gu+Aout*(fi).^2.*(u-uout)-sifiu*lapfi;
       
       Fu=-Av*duu*lapu+Gu;
   
        H=Fu;
        lap3D
        lapFu=lapH;
        
         Ft=-sifiu*lapu; %  surface tension between membrane and u
%         H=Ft;
%         lap3D
%         lapFt=lapH;        
%         
%         %   dynamical equations,  for conservation of mass use  fi=fi-dt*(F+Fs);
         I=20*Fu*sum(sum(sum((fi>=-.99))))/Nx/Ny/Nz;
         I(fi<=-.99)=0;
         I;
         H=F-sigma*Fs+Ft;
         lap3D
         lapE=lapH;
%fi=fi-dt*(F+Fs+Ft);
        fi=fi+dt*(lapE+I);
        fi(:,:,1)=fi(:,:,2);
       % Iu=sum(sum(sum(u)))/Nx/Ny/Nz-u0;
        u=u+Du*dt*(lapFu);
    end
    ux(:,:)=u(:,Ny/2,:);
    
    %%   
h=max(max(max(isnan(fi(:,:,:)))));
    if h==1;
        h
        break
    end
      
    iter
    
    %%
                
figure(1)
clf
%fim = smooth3(fi,'box',5);
isosurface(fi,0)
%camlight(10,40,'infinity')
view(64,26),
axis equal
axis([1 Nx 1 Ny 1 Nz])
light
material metal
%%
    figure(2)

%     u = smooth3(u,'box',3); 
    clf
       
    %cdata = smooth3((u-min(min(min(u))))./(max(max(max(u)))-min(min(min(u)))),'box',5);
       cdata = smooth3(u,'box',5);
    fim = smooth3(fi,'box',3);
    p4=patch(isosurface(fim,0));
    isonormals(fim,p4);
    isocolors(cdata,p4);
    set(p4,'FaceColor','interp','EdgeColor','none'),
    camlight, lighting phong
    axis equal, view(-14,40), axis off
    axis([1 Nx 1 Ny 1 Nz/2]),
    light
material metal
    colorbar

% cdata = smooth3((u-min(min(min(u))))./(max(max(max(u)))-min(min(min(u)))),'box',5);
% [x,y,z] = meshgrid(1:1:Nx,1:1:Ny,1:1:Nz);
% xslice = [Nx/4:Nx/4:3*Nx/4];yslice = [Ny/4:Ny/4:3*Ny/4]; zslice = [R:R:3*R];
% p3=slice(x,y,z,u,xslice,yslice,zslice);
% set(p3,'FaceColor','interp','EdgeColor','none','FaceAlpha',0.5);
% axis equal, view(73,8), 
% colorbar

%para hacer peliculas 3D usar cine3D y cine2D

    Fm(:,:,:,iter)=fi(:,:,:);
    U(:,:,:,iter)=u(:,:,:);
    %%
    figure(3)
    fix(:,:)=fi(:,Ny/2,:);
    contour(fix,[0 0],'k')
    hold on
    contour(fix0,[0 0],'r')
    axis equal
    %getframe(gcf);
    hold off
    
end

