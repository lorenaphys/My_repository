% Program to calculate phase fiels in 3 dimensions

clear all

et=.01;
dx=1;
NF=800;
sig=0*(1:NF);
ep1=2;
ep=ep1^2;
sigma=-.1;
Nx=40;
Ny=40;
Nz=60;
R=9;
N=2;
sifiu=0.;
duu=.1; 
Du=10;
Dfi=1;
eta=5;


%%%%%%%%%%%%%%%%%%%%% strength of the fields  %%%%%%%%%%%
Aout=2; 
Av=2; 
Afi=.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% form of fi and initial values

semiesf3D
%prisma3D
fiini=fi;
%%%%%%%%%%% parameters for iteraion loop %%%%%%%%%%%%%%%%%%%%%%%%
step=200;
iter=1;
dt=1e-5;
cont=iter;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initial conditions for u and beta %%%%%

% u1=1; 
% u2=-.0; 
% 
% uout=0.0; 

bet=.1;

% uinit=0*fi;
% [uf vf]=find(fi>=-.9);
% for i=1:length(uf)
%  uinit(uf(i),vf(i))=u0+.1*(rand-.5);%+u(uf(i),vf(i));
% end
[X,Y,Z]=meshgrid(1:Nx,1:Ny,1:Nz);
        teta=atan2((Y-Ny/2),(X-Nx/2));
        rad=sqrt((X-Nx/2+.5).^2+(Y-Ny/2+.5).^2);
        %u=1.5*exp(-((X-Nx/2-.5).^2+(Y-Ny/2-.5).^2+(Z-R+2).^2)/20);
        u=2.5*rad.*(cos(teta*N)+sin(teta*N)).*(Z/Nz)/max(max(max(rad)))+(exp(-((-X+Nx/2).^2+(-Y+Ny/2).^2+(-Z+(R+14)).^2)/80));
        %u=2.5*rad.*(cos(teta*N)+sin(teta*N)).*(Z/Nz)/max(max(max(rad)))+(exp(-((X-Nx/3).^2+(Y-Ny/3).^2+(Z-(R+2)).^2)/50));
        %u=2.5*rand(Nx,Ny,Nz);
%u=1-u;
        %u(fi<=-.9)=0;
        u0=sum(sum(sum(u)))/Nx/Ny/Nz;
 %       u(fi>=.9)=0;        
%u=smooth3(uinit,'gaussian',3);



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
%%    store the initial domain

%fi = smooth3(fi,'box',5);
fix0(:,:)=fi(Nx/2,:,:);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

%%
%load phase3d-2
% clear all
% load isoformas-may27-1
% NF=800;
% cont=iter;
for iter=cont:NF              %time loop
    for iiter=1:step
 
 %%   delficiones 
 
        H=fi;
        lap3D
        lapfi=lapH;
        
        mu=((fi-ep1.*(bet*u.^2)).*((fi).^2-1)-ep*lapfi);        
        
        H=mu;
        lap3D
        lapmu=lapH;
        
        F=Afi*((3*fi.^2-1-2*ep1*fi.*(bet*u.^2)).*mu-ep*lapmu);     
        
        
        H=F;
        lap3D
        lapF=lapH;

        H=u;
        lap3D
        lapu=lapH;
        
        Fs=lapfi;
        
        
        H=Fs;
        lap3D
        lapFs=lapH;


%%   ahora  la u  %%%%%%%%%%%%%


       Gu=-Afi*ep1*bet*u.*mu.*((fi).^2-1)-sifiu*lapfi; 

       
       Fu=-Av*duu*lapu+Gu;
       gFu=grad3DR(Fu);
        H=Fu;
        lap3D
        lapFu=lapH;
        
        Ft=-sifiu*lapu; %  surface tension between membrane and u
        H=Ft;
        lap3D
        lapFt=lapH; 
%%   tensor de esfurezos
        gfi=grad3DR(fi);
        gmu=grad3DR(mu);
        P=(Afi*mu.^2-sigma*abs(gfi(:,:,:,1).^2+gfi(:,:,:,2).^2+gfi(:,:,:,3).^2)-fi(:,:,:).*F(:,:,:));
        ggfi1=grad3DR(gfi(:,:,:,1));
        ggfi2=grad3DR(gfi(:,:,:,2));
        ggfi3=grad3DR(gfi(:,:,:,3));        
        

str=zeros(Nx,Ny,Nz,3,3);
 str(:,:,:,1,1)=eta*(P(:,:,:)+2*sigma*gfi(:,:,:,1).*gfi(:,:,:,1)-2*Afi*ep*gmu(:,:,:,1).*gfi(:,:,:,1)-2*Afi*ep*ggfi1(:,:,:,1));
 str(:,:,:,2,2)=eta*(P(:,:,:)+2*sigma*gfi(:,:,:,2).*gfi(:,:,:,2)-2*Afi*ep*gmu(:,:,:,2).*gfi(:,:,:,2)-2*Afi*ep*ggfi2(:,:,:,2));
 str(:,:,:,3,3)=eta*(P(:,:,:)+2*sigma*gfi(:,:,:,3).*gfi(:,:,:,3)-2*Afi*ep*gmu(:,:,:,3).*gfi(:,:,:,3)-2*Afi*ep*ggfi3(:,:,:,3));
 
 str(:,:,:,1,2)=eta*(2*sigma*gfi(:,:,:,1).*gfi(:,:,:,2)-2*Afi*ep*gmu(:,:,:,1).*gfi(:,:,:,2)-2*Afi*ep*ggfi2(:,:,:,1));
 str(:,:,:,1,3)=eta*(2*sigma*gfi(:,:,:,1).*gfi(:,:,:,3)-2*Afi*ep*gmu(:,:,:,1).*gfi(:,:,:,3)-2*Afi*ep*ggfi3(:,:,:,1));
 str(:,:,:,2,1)=eta*(2*sigma*gfi(:,:,:,2).*gfi(:,:,:,1)-2*Afi*ep*gmu(:,:,:,2).*gfi(:,:,:,2)-2*Afi*ep*ggfi1(:,:,:,2));
 str(:,:,:,2,3)=eta*(2*sigma*gfi(:,:,:,2).*gfi(:,:,:,3)-2*Afi*ep*gmu(:,:,:,2).*gfi(:,:,:,3)-2*Afi*ep*ggfi3(:,:,:,2));
 str(:,:,:,3,1)=eta*(2*sigma*gfi(:,:,:,3).*gfi(:,:,:,1)-2*Afi*ep*gmu(:,:,:,3).*gfi(:,:,:,1)-2*Afi*ep*ggfi1(:,:,:,3));
 str(:,:,:,3,2)=eta*(2*sigma*gfi(:,:,:,3).*gfi(:,:,:,2)-2*Afi*ep*gmu(:,:,:,3).*gfi(:,:,:,2)-2*Afi*ep*ggfi2(:,:,:,3)); 

 
 
 dS1(:,:,:)=str(:,:,:,1,1).*gFu(:,:,:,1)+str(:,:,:,1,2).*gFu(:,:,:,2)+str(:,:,:,1,3).*gFu(:,:,:,3);
 dS2(:,:,:)=str(:,:,:,2,1).*gFu(:,:,:,1)+str(:,:,:,2,2).*gFu(:,:,:,2)+str(:,:,:,2,3).*gFu(:,:,:,3);
 dS3(:,:,:)=str(:,:,:,3,1).*gFu(:,:,:,1)+str(:,:,:,3,2).*gFu(:,:,:,2)+str(:,:,:,3,3).*gFu(:,:,:,3);
 
 gs1=grad3DR(dS1(:,:,:));
 gs2=grad3DR(dS2(:,:,:));
 gs3=grad3DR(dS3(:,:,:));
 S=gs1(:,:,:,1)+gs2(:,:,:,2)+gs3(:,:,:,3);
 
%%            dynamical equations,  for conservation of mass use  fi=fi-dt*(F+Fs);

         I=200.*u*sum(sum(sum((fi>=-.99))))/Nx/Ny/Nz;
         I(find(abs(fi)>=.9))=0;
         I;
         H=F-sigma*Fs+Ft;
         lap3D
         lapE=lapH;
%fi=fi-dt*(F+Fs+Ft);
        fi=fi+Dfi*dt*(lapE+I);
        fi(:,:,1)=fi(:,:,2);
       % Iu=sum(sum(sum(u)))/Nx/Ny/Nz-u0;
        u=u+dt*(Du*lapFu+S);
        u(:,:,1)=u(:,:,2);        
    end
    ux(:,:)=u(:,Ny/2,:);
    fix(:,:)=fi(Nx/2,:,:);

    %%   
h=max(max(max(isnan(fi(:,:,:)))));
    if h==1;
       'nans'
        break
    end
      
    iter
    
    %%
                
figure(1)
clf
%fim = smooth3(fi,'box',5);
%hold on
%isosurface(S)
mesh(fix)
%camlight(10,40,'infinity')
 view(25,40),
% axis equal
% axis([1 Nx 1 Ny 1 Nz])
light
%material metal
hold off
%%
      figure(2)
     clf
[x,y,z] = meshgrid(1:1:Ny,1:1:Nx,1:1:Nz);
xslice = [Nx/2-R:R:Nx/2+R,Nx/2-R:R:Nx/2+R];yslice = [Ny/2-R:R:Ny/2+R,Ny/2-R:R:Ny/2+R]; zslice = [0:3:R,0:3:R];
p3=slice(x,y,z,S,xslice,yslice,zslice);
set(p3,'FaceColor','flat','EdgeColor','none','FaceAlpha',0.1);
rs=max(abs(max(max(max(S)))),abs(min(min(min(S)))));

axis equal, view(74,18), 
set(gca,'CLim',[-rs,rs])
colorbar;
%%
    figure(3)
    
    

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
    axis([1 Nx 1 Ny 1 Nz]),
    light
material metal
    colorbar



%para hacer peliculas 3D usar cine3D y cine2D
%%
    Fm(:,:,:,iter)=fi(:,:,:);
    U(:,:,:,iter)=u(:,:,:);
    Sm(:,:,:,iter)=S(:,:,:);
    %%
    figure(4)
    fix(:,:)=fi(:,Ny/2,:);
    contour(fix,[0 0],'k')
    hold on
    contour(fix0,[0 0],'r')
    axis equal
    %getframe(gcf);
    hold off
    
end

