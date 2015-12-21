%% Program to calculate phase fiels and Turing-Stress interactions in 3 dimensions

clear all


%%  form of fi and initial values
NF=400;
Nx=30;
Ny=30;
Nz=30;
R=9;
N=0;

semiesf3D
%prisma3D
%% %%%%%%%%% parameters for iteraion loop %%%%%%%%%%%%%%%%%%%%%%%%
step=200;
iter=1;
dt=1e-5;
cont=iter;


%% %%%%%%% inital parameter for Turing 
etaT=sqrt(2);
Dus=1;
Du=Dus*.516/etaT;
Dv=Dus/etaT;
%Du=.2/etaT;

hT=-1.;
aT=1/.899;
bT=-.91/.899;
cT=0.02;

dtT=500*dt;

%% %%%%%%%%%%%%%%%%%%% strength of the fields  %%%%%%%%%%%
ep1=2;
ep=ep1^2;
sigma=-.01;
bet=.2;
sifiu=0.;
duu=.1; 
%Du=10;
Dfi=1;
eta=0.01;%Du/2;
As=2; 
Av=2; 
Afi=.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initial conditions for u,v and beta %%%%%

[X,Y,Z]=meshgrid(1:Nx,1:Ny,1:Nz);
        
        tetar=0;   % X rotation
        fir=0;     % Z rotation
        RX=(X)*cos(fir)-(Y)*sin(fir)*cos(tetar)+(Z)*sin(fir)*sin(tetar);
        RY=(X)*sin(fir)+(Y)*cos(fir)*cos(tetar)-(Z)*cos(fir)*sin(tetar);
        RZ=(Y)*sin(tetar)+(Z)*cos(tetar);
        teta=atan2((RY-Ny/2),(RX-Nx/2));
        rad=sqrt((RX-Nx/2+.25).^2+(RY-Ny/2+.25).^2);
        u=2.5*rad.*(cos(teta*N)+sin(teta*N)).*(RZ/Nz)/max(max(max(rad)))+(exp(-((RX-Nx/2).^2+(RY-Ny/2).^2+(RZ-(R+2)).^2)/20));
        u0=sum(sum(sum(u)))/Nx/Ny/Nz;
        
       
% 
%     figure(3)
%     clf
%        
%     %cdata = smooth3((u-min(min(min(u))))./(max(max(max(u)))-min(min(min(u)))),'box',5);
%        cdata = smooth3(u,'box',5);
%     fim = smooth3(fi,'box',3);
%     p4=patch(isosurface(fim,0));
%     isonormals(fim,p4);
%     
%     isocolors(cdata,p4);
%     set(p4,'FaceColor','interp','EdgeColor','none'),
%     camlight, lighting phong
%     axis equal, view(-14,40), axis off
%     axis([1 Nx 1 Ny 1 Nz/2]),
%     light
%     material metal
%     colorbar
%     
%           figure(2)
%      clf
% [x,y,z] = meshgrid(1:1:Ny,1:1:Nx,1:1:Nz);
% xslice = [Nx/2-R:R:Nx/2+R,Nx/2-R:R:Nx/2+R];yslice = [Ny/2-R:R:Ny/2+R,Ny/2-R:R:Ny/2+R]; zslice = [0:3:R,0:3:R];
% p3=slice(x,y,z,u,xslice,yslice,zslice);
% set(p3,'FaceColor','flat','EdgeColor','none','FaceAlpha',0.1);
% rs=max(abs(max(max(max(u)))),abs(min(min(min(u)))));
% 
% axis equal, view(74,18), 
% set(gca,'CLim',[-rs,rs])
% colorbar;
%%    store the initial domain

    v=0.1*(rand(Nx,Ny,Nz)-.5);        

    u(fi<=-.9)=0;
    v(fi<=-.9)=0;
    fiini=fi;
    fix0(:,:)=fi(Nx/2,:,:);

%%
%load phase3d-2
% clear all
% load isoformas-may27-1
% NF=800;
% cont=iter;
%% DYNAMICS AND TURING
for iter=cont:NF              %time loop
    for iiter=1:step
 
 %%   Phase field fi
 
        H=fi;
        lap3D
        lapfi=lapH;
        
        mu=((fi-ep1.*(bet*u.^2)).*((fi).^2-1)-ep*lapfi);        
        
        H=mu;
        lap3D
        lapmu=lapH;

        H=u;
        lap3D
        lapu=lapH;
        
        F=Afi*2*((3*fi.^2-1-2*ep1*fi.*(bet*u.^2)).*mu-ep*lapmu)-sifiu*lapu;     
         
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


%%  Phase field  u


       Gu=-Afi*4*ep1*bet*u.*mu.*((fi).^2-1)-sifiu*lapfi; 

       
       Fu=-Av*duu*lapu+Gu;
       gFu=grad3DR(Fu);
        H=Fu;
        lap3D
        lapFu=lapH;
        
        Ft=-sifiu*lapu; %  surface tension between membrane and u
        H=Ft;
        lap3D
        lapFt=lapH;
        
%%  Stress tensor
        gfi=grad3DR(fi);
        gmu=grad3DR(mu);
        gu=grad3DR(u);
        P=(Afi*mu.^2-sigma*abs(gfi(:,:,:,1).^2+gfi(:,:,:,2).^2+gfi(:,:,:,3).^2)+0.5*Av*duu*abs(gu(:,:,:,1).^2+gu(:,:,:,2).^2+gu(:,:,:,3).^2)+0.5*As*u.^2+sifiu*(gfi(:,:,:,1).*gu(:,:,:,1)+gfi(:,:,:,2).*gu(:,:,:,2)+gfi(:,:,:,3).*gu(:,:,:,3))-fi(:,:,:).*(F(:,:,:)-sigma*Fs(:,:,:)+Ft(:,:,:)));
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
 
%%  %%%%%%%%%%%  Dynamical equations,  for conservation of mass use  fi=fi-dt*(F+Fs);

         I=200.*u*sum(sum(sum((fi>=-.99))))/Nx/Ny/Nz;
         I(find(abs(fi)>=.9))=0;
         I;
         H=F-2*sigma*Fs+Ft;
         lap3D
         lapE=lapH;
        %fi=fi-dt*(F+Fs+Ft);
        fi=fi+Dfi*dt*(lapE+I);
        fi(:,:,1)=fi(:,:,2);
        % Iu=sum(sum(sum(u)))/Nx/Ny/Nz-u0;
        %u=u+dt*(Du*lapFu+S);
% %%%%%%%%%  Turing interaction Laplacoano con |nabla^2=oen la frontera del cubo 
		for tir=1:10
        
		H=u;
        lap3Dt
        lapu=lapH;
        
        H=v;
        lap3Dt
        lapv=lapH; 
        			
        lapu(fi<=-.9)=0;
        lapv(fi<=-.9)=0;
        u(fi<=-.99)=0;
        v(fi<=-.99)=0;

        u=u+dtT*(Du*(As*lapu+lapFu)+S+(u+aT*v-cT*u.*v-u.*v.^2));
        v=v+dtT*(Dv*lapv+(bT*v+hT*u+cT*u.*v+u.*v.^2));        
        
		end 
      
        %u(:,:,1)=u(:,:,2); 
        %v(:,:,1)=v(:,:,2);
    end
    ux(:,:)=u(:,Ny/2,:);
    vx(:,:)=v(:,Ny/2,:);
    fix(:,:)=fi(Nx/2,:,:);

    %%   
h=max(max(max(isnan(fi(:,:,:)))));
    if h==1;
       'nans'
        break
    end
      
    iter

%%  For movies in 3 dimension
    Fim(:,:,:,iter)=fi(:,:,:);
    Um(:,:,:,iter)=u(:,:,:);
    Vm(:,:,:,iter)=v(:,:,:);
    Sm(:,:,:,iter)=S(:,:,:);
    %%
                
figure(1)
clf
%fim = smooth3(fi,'box',5);
%hold on
%isosurface(S)
%mesh(fix)
hold on
surf(fix,ux,'FaceAlpha',0.8),shading interp, grid
%camlight(10,40,'infinity')
 view(25,40),
% axis equal
% axis([1 Nx 1 Ny 1 Nz])
light
%material metal
hold off
%%
%       figure(2)
%      clf
% [x,y,z] = meshgrid(1:1:Ny,1:1:Nx,1:1:Nz);
% xslice = [Nx/2-R:R:Nx/2+R,Nx/2-R:R:Nx/2+R];yslice = [Ny/2-R:R:Ny/2+R,Ny/2-R:R:Ny/2+R]; zslice = [0:3:R,0:3:R];
% p3=slice(x,y,z,S,xslice,yslice,zslice);
% set(p3,'FaceColor','flat','EdgeColor','none','FaceAlpha',0.1);
% rs=max(abs(max(max(max(S)))),abs(min(min(min(S)))));
% 
% axis equal, view(74,18), 
% set(gca,'CLim',[-rs,rs])
% colorbar;

%%
%      figure(3)
%      clf
% p3=slice(x,y,z,u,xslice,yslice,zslice);
% set(p3,'FaceColor','flat','EdgeColor','none','FaceAlpha',0.1);
% rsu=max(abs(max(max(max(u)))),abs(min(min(min(u)))));
% 
% axis equal, view(74,18), 
% set(gca,'CLim',[-rsu,rsu])
% colorbar;
%%
    figure(4)
    
    clf

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

%%
    figure(5)
    fix(:,:)=fi(:,Ny/2,:);
    contour(fix,[0 0],'k')
    hold on
    contour(fix0,[0 0],'r')
    axis equal
    %getframe(gcf);
    hold off
    
%%    
    
%     figure(6)
% clf
% hold on
% surf(fix,vx,'FaceAlpha',0.8),shading interp, grid
% view(25,40),
% light    
end

