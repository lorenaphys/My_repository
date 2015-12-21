% Program to calculate phase fiels in 3 dimensions

clear all


dx=1;
NF=200;
sig=0*(1:NF);
ep1=2;
ep=ep1^2;
sigma=-.1;
Nx=30;
Ny=30;
Nz=30;
R=11;
N=6;
sifiu=0.;
duu=.1; 
Dus=1;
Dfi=1;
eta=.01;

%% %%%% Parametros del Turing
 eta1=sqrt(3);
 Du=Dus*.516/eta1;
 Dv=Dus/eta1;
h=-1.;
a=1/.899;
b=-.91/.899;
c=0.02;
dt1=.02;
cT=1;  %tiempo en el que empieza el turing
%%
%%%%%%%%%%%%%%%%%%%%% strength of the fields  %%%%%%%%%%%
Aout=2; 
Av=2; 
Afi=.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% form of fi and initial values
%load dominio
%semiesf3D
formas3D
%prism3D
fiini=fi;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initial conditions for u and beta %%%%%

% u1=1; 
% u2=-.0; 
% 
% uout=0.0; 

bet=.1;

        
[X,Y,Z]=meshgrid(1:Nx,1:Ny,1:Nz);
%         
%         tetar=0;   % X rotation
%         fir=0;     % Z rotation
%         RX=(X)*cos(fir)-(Y)*sin(fir)*cos(tetar)+(Z)*sin(fir)*sin(tetar);
%         RY=(X)*sin(fir)+(Y)*cos(fir)*cos(tetar)-(Z)*cos(fir)*sin(tetar);
%         RZ=(Y)*sin(tetar)+(Z)*cos(tetar);
%         teta=atan2((RY-Ny/2),(RX-Nx/2));
%         rad=sqrt((RX-Nx/2+.25).^2+(RY-Ny/2+.25).^2);
%        u=2.5*rad.*(cos(teta*N)+sin(teta*N)).*(RZ/Nz)/max(max(max(rad)))+(exp(-((RX-Nx/2).^2+(RY-Ny/2).^2+(RZ-(R+2)).^2)/20));
        u=.05*(1.5*(exp(-((X-Nx/2).^2+(Y-Ny/2).^2+(Z-(R+2)).^2)/10))-.1);
  u=u+.2*(rand(Nx,Ny,Nz)-.5);
        v=u+.2*(rand(Nx,Ny,Nz)-.5);       
% 
%         u(fi<=-.9)=0;
 %        v(fi<=-.9)=0;
%%    store the initial domain

%fi = smooth3(fi,'box',5);
%load temp
fix0(:,:)=fi(Nx/2,:,:);
%%%%%%%%%%% parameters for iteraion loop %%%%%%%%%%%%%%%%%%%%%%%%
step=20;
iter=1;
dt=1e-5;
cont=iter;
ct=0;
%%

for iter=cont:NF  
    iter%time loop
    for iiter=1:step
 
 %%   delficiones 
 
        H=fi;
        lap3Dt
        lapfi=lapH;
        
        mu=((fi-ep1.*(bet*u.^2)).*((fi).^2-1)-ep*lapfi);        
        
        H=mu;
        lap3Dt
        lapmu=lapH;
        
        F=Afi*((3*fi.^2-1-2*ep1*fi.*(bet*u.^2)).*mu-ep*lapmu);     
        
        
        H=F;
        lap3Dt
        lapF=lapH;

        H=u;
        lap3Dt
        lapu=lapH;
        
        Fs=lapfi;
        
        
        H=Fs;
        lap3Dt
        lapFs=lapH;


%%   ahora  la u  %%%%%%%%%%%%%
% if iter == cT
%          u=.1*u+.2*(rand(Nx,Ny,Nz)-.5);
%         v=.1*u+.2*(rand(Nx,Ny,Nz)-.5);
% end
% if iter<=cT

            uts
%        end
%%        


       Gu=-Afi*ep1*bet*u.*mu.*((fi).^2-1)-sifiu*lapfi; 

       
       Fu=-Av*duu*lapu+Gu;
      gFu=grad3DR(Fu);
        H=Fu;
        lap3Dt
        lapFu=lapH;
        
         Ft=-sifiu*lapu; %  surface tension between membrane and u
        H=Ft;
        lap3Dt
        lapFt=lapH; 
%%   tensor de esfuerzos
        gfi=grad3DR(fi);
        gmu=grad3DR(mu);
        gu=grad3DR(u);
        P=(Afi*mu.^2-sigma*abs(gfi(:,:,:,1).^2+gfi(:,:,:,2).^2+gfi(:,:,:,3).^2)+0.5*Av*duu*abs(gu(:,:,:,1).^2+gu(:,:,:,2).^2+gu(:,:,:,3).^2)+sifiu*(gfi(:,:,:,1).*gu(:,:,:,1)+gfi(:,:,:,2).*gu(:,:,:,2)+gfi(:,:,:,3).*gu(:,:,:,3))-fi(:,:,:).*(F(:,:,:)-sigma*Fs(:,:,:)+Ft(:,:,:)));
        ggfi1=grad3DR(gfi(:,:,:,1));
        ggfi2=grad3DR(gfi(:,:,:,2));
        ggfi3=grad3DR(gfi(:,:,:,3));        
        

    str=zeros(Nx,Ny,Nz,3,3);
    str(:,:,:,1,1)=eta*(P(:,:,:)+2*sigma*gfi(:,:,:,1).*gfi(:,:,:,1)-Av*duu*gfi(:,:,:,1).*gu(:,:,:,1)-2*Afi*ep*gmu(:,:,:,1).*gfi(:,:,:,1)+2*Afi*ep*mu.*ggfi1(:,:,:,1));
    str(:,:,:,2,2)=eta*(P(:,:,:)+2*sigma*gfi(:,:,:,2).*gfi(:,:,:,2)-Av*duu*gfi(:,:,:,2).*gu(:,:,:,2)-2*Afi*ep*gmu(:,:,:,2).*gfi(:,:,:,2)+2*Afi*ep*mu.*ggfi2(:,:,:,2));
    str(:,:,:,3,3)=eta*(P(:,:,:)+2*sigma*gfi(:,:,:,3).*gfi(:,:,:,3)-Av*duu*gfi(:,:,:,3).*gu(:,:,:,3)-2*Afi*ep*gmu(:,:,:,3).*gfi(:,:,:,3)+2*Afi*ep*mu.*ggfi3(:,:,:,3));
 
    str(:,:,:,1,2)=eta*(2*sigma*gfi(:,:,:,1).*gfi(:,:,:,2)-Av*duu*gfi(:,:,:,1).*gu(:,:,:,2)-2*Afi*ep*gmu(:,:,:,1).*gfi(:,:,:,2)+2*Afi*ep*mu.*ggfi1(:,:,:,2));
    str(:,:,:,1,3)=eta*(2*sigma*gfi(:,:,:,1).*gfi(:,:,:,3)-Av*duu*gfi(:,:,:,1).*gu(:,:,:,3)-2*Afi*ep*gmu(:,:,:,1).*gfi(:,:,:,3)+2*Afi*ep*mu.*ggfi1(:,:,:,3));
    str(:,:,:,2,1)=eta*(2*sigma*gfi(:,:,:,2).*gfi(:,:,:,1)-Av*duu*gfi(:,:,:,2).*gu(:,:,:,1)-2*Afi*ep*gmu(:,:,:,2).*gfi(:,:,:,1)+2*Afi*ep*mu.*ggfi2(:,:,:,1));
    str(:,:,:,2,3)=eta*(2*sigma*gfi(:,:,:,2).*gfi(:,:,:,3)-Av*duu*gfi(:,:,:,2).*gu(:,:,:,3)-2*Afi*ep*gmu(:,:,:,2).*gfi(:,:,:,3)+2*Afi*ep*mu.*ggfi2(:,:,:,3));
    str(:,:,:,3,1)=eta*(2*sigma*gfi(:,:,:,3).*gfi(:,:,:,1)-Av*duu*gfi(:,:,:,3).*gu(:,:,:,1)-2*Afi*ep*gmu(:,:,:,3).*gfi(:,:,:,1)+2*Afi*ep*mu.*ggfi3(:,:,:,1));
    str(:,:,:,3,2)=eta*(2*sigma*gfi(:,:,:,3).*gfi(:,:,:,2)-Av*duu*gfi(:,:,:,3).*gu(:,:,:,2)-2*Afi*ep*gmu(:,:,:,3).*gfi(:,:,:,2)+2*Afi*ep*mu.*ggfi3(:,:,:,2)); 
        
 
 
 dS1(:,:,:)=str(:,:,:,1,1).*gFu(:,:,:,1)+str(:,:,:,1,2).*gFu(:,:,:,2)+str(:,:,:,1,3).*gFu(:,:,:,3);
 dS2(:,:,:)=str(:,:,:,2,1).*gFu(:,:,:,1)+str(:,:,:,2,2).*gFu(:,:,:,2)+str(:,:,:,2,3).*gFu(:,:,:,3);
 dS3(:,:,:)=str(:,:,:,3,1).*gFu(:,:,:,1)+str(:,:,:,3,2).*gFu(:,:,:,2)+str(:,:,:,3,3).*gFu(:,:,:,3);
 
 gs1=grad3DR(dS1(:,:,:));
 gs2=grad3DR(dS2(:,:,:));
 gs3=grad3DR(dS3(:,:,:));
 S=gs1(:,:,:,1)+gs2(:,:,:,2)+gs3(:,:,:,3);
 
%%            dynamical equations,  for conservation of mass use  fi=fi-dt*(F+Fs);

         I=200.*u;
         I(find(abs(fi)>=.9))=0;
         I;
         H=F-2*sigma*Fs+Ft;
         lap3Dt
         lapE=lapH;
%fi=fi-dt*(F+Fs+Ft);
        fi=fi+Dfi*dt*(lapE+I);
        fi(:,:,1)=fi(:,:,2);
       % Iu=sum(sum(sum(u)))/Nx/Ny/Nz-u0;
        u=u+dt*(Dus*lapFu+S);
        %uts
       % u(:,:,1)=u(:,:,2);  
           
    end
    ux(:,:)=u(:,Ny/2,:);
    fix(:,:)=fi(Nx/2,:,:);

    %%   
hh=max(max(max(isnan(fi(:,:,:)))));
    if hh==1;
       'nans'
        break
    end
      

    
    %%
                
figure(1)
clf
    ux(:,:)=u(:,Ny/2,:);
    vx(:,:)=v(:,Ny/2,:);
    fix(:,:)=fi(Nx/2,:,:);
clf
hold on
mesh(ux,'FaceAlpha',0.5),shading interp,  view(-36,18)
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
    figure(4)
    fix(:,:)=fi(:,Ny/2,:);
    contour(fix,[0 0],'k')
    hold on
    contour(fix0,[0 0],'r')
    axis equal
    %getframe(gcf);
    hold off
    
end

