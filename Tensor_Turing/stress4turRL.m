% Program to calculate phase fiels in 3 dimensions

%clear all


%dx=1;
NF=100;
%sig=0*(1:NF);
ep1=2;
ep=ep1^2;
sigma=-0.1;
Nx=40;
Ny=40;
Nz=70;
sifiu=0.;
%duu=.1; 
Du=20;
Dfi=1;
eta=1;
Afi = 0.5;
%Ab = 0.5;
Af = 0.1;
As = 0.1;
beta = 0.1;
u1 = 0;
u2 = 1;
u3 = 0;
L = 1;

%% %%%% Parametros del Turing
 %eta1=sqrt(3);
 %Du=Dus*.516/eta1;
 %Dv=Dus/eta1;
%h=-1.;
%a=1/.899;
%b=-.91/.899;
%c=0.02;
%dt1=.02;
%cT=1;  %tiempo en el que empieza el turing
%%
%%%%%%%%%%%%%%%%%%%%% strength of the fields  %%%%%%%%%%%
%Aout=2; 
%Av=2; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% form of fi and initial values
%load dominio

% load initial
% u0=uout+ues; v0=vout+ves; w0=wout+wes; x0=xout+xes;
fi=ones(Nx,Ny,Nz);
r = zeros(Nx,Ny,Nz);

for i=1:Nx
    for j=1:Ny
        for k=1:Nz
      r(i,j,k)=sqrt((i-Nx/2)^2+(j-Ny/2)^2+(k)^2);
      if r(i,j,k)>=15
      fi(i,j,k)=-1;
      end
        end
   
   end
end
%formas3D
%prism3D
%fiini=fi;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initial conditions for u and beta %%%%%

% u1=1; 
% u2=-.0; 
% 
% uout=0.0; 

%bet=.1;

    %R = 9;
    %N = 6;
    %[X,Y,Z]=meshgrid(1:Nx,1:Ny,1:Nz);
        %teta=atan2((Y-Ny/2),(X-Nx/2));
        %rad=sqrt((X-Nx/2+.5).^2+(Y-Ny/2+.5).^2);
        %u=1.5*exp(-((X-Nx/2-.5).^2+(Y-Ny/2-.5).^2+(Z-R+2).^2)/100);
        %u=2.5*rad.*(cos(teta*N)+sin(teta*N)).*(Z/Nz)/max(max(max(rad)))+(exp(-((-X+Nx/2-.5).^2+(-Y+Ny/2-.5).^2+(-Z+R+14).^2)/80));
         %u=-2.5*rad.*(cos(teta*N)+sin(teta*N)).*(Z/Nz)/max(max(max(rad)))+(exp(-((X-Nx/2).^2+(Y-Ny/2).^2+(Z-R+14).^2)/50));
        %u=2.5*rand(Nx,Ny,Nz);
    
u = zeros(Nx,Ny,Nz);
        
     for i=1:Nx
         for j=1:Ny
             for k=1:Nz
                 u(i,j,k)=1.5*exp(-((i-Nx/2-0.5)^2+(j-Ny/2-0.5)^2+(k-7-0.5)^2)/50);%primera prueba
             end
         end
     end        

Fm = zeros(Nx,Ny,Nz,NF+1);
Um = zeros(Nx,Ny,Nz,NF+1);
%Sm = zeros(Nx,Ny,Nz,NF);

Fm(:,:,:,1) = fi;
Um(:,:,:,1) = u;

%fix0(:,:)=fi(Nx/2,:,:);
%%%%%%%%%%% parameters for iteraion loop %%%%%%%%%%%%%%%%%%%%%%%%
step=10;
dt=1e-5;
%ct=0;
%%

t = tic();

for iter = 1:NF  
    %iter%time loop
    for iiter = 1:step
 
 %%   deficiones 
 
		lapfi = lapf3D(fi);
        
        mu=(fi-ep1*beta*u).*((fi).^2-1)-ep*lapfi;        
        
		lapmu = lapf3D(mu);
        
        F=Afi*((3*fi.^2-1-2*ep1*beta*fi.*u).*mu-2*ep*lapmu);     
        
		lapmu = lapf3D(mu);

		lapu = lapf3D(u);
        
        Fs=lapfi;
        
        %lapFs = lapf3D(Fs);


%%   ahora  la u  %%%%%%%%%%%%%
% if iter == cT
%          u=.1*u+.2*(rand(Nx,Ny,Nz)-.5);
%         v=.1*u+.2*(rand(Nx,Ny,Nz)-.5);
% end
% if iter<=cT

            %uts
%        end
%%        

	Gup = 2*As*(fi.^2-1).*(u-u1).*(u-u2).*(2*u-u1-u2)+2*Af*fi.^2.*(u-u3);
        Gu=-2*Afi*ep1*beta*mu.*((fi).^2-1)-sifiu*lapfi+Gup; 

       
        Fu=-2*As*L*lapu+Gup+Gu;
        gFu=grad3DR(Fu);
       
	lapFu = lapf3D(Fu);
        
        Ft=-sifiu*lapu; %  surface tension between membrane and u
        
	lapFt = lapf3D(Ft);
        
        F1 = 4*As*fi.*(fi.^2-1).*(u-u1).^2.*(u-u2).^2 +2*Af*fi.*(u-u3).^2;
        
	lapF1 = lapf3D(F1);
        
%%   tensor de esfuerzos
        gfi=grad3DR(fi);
        gmu=grad3DR(mu);
        gu=grad3DR(u);
        P=(Afi*mu.^2-sigma*abs(gfi(:,:,:,1).^2+gfi(:,:,:,2).^2+gfi(:,:,:,3).^2)+abs(gu(:,:,:,1).^2+gu(:,:,:,2).^2+gu(:,:,:,3).^2)...
           +sifiu*(gfi(:,:,:,1).*gu(:,:,:,1)+gfi(:,:,:,2).*gu(:,:,:,2)+gfi(:,:,:,3).*gu(:,:,:,3)) +As*(fi.^2-1).^2.*(u-u1).^2.*(u-u2).^2 ...
           +Af*fi.^2.*(u-u3).^2-fi(:,:,:).*(F(:,:,:)-2*sigma*Fs(:,:,:)+Ft(:,:,:)+F1(:,:,:)));
        ggfi1=grad3DR(gfi(:,:,:,1));
        ggfi2=grad3DR(gfi(:,:,:,2));
        ggfi3=grad3DR(gfi(:,:,:,3));        
        

    str=zeros(Nx,Ny,Nz,3,3);
    str(:,:,:,1,1)=eta*(P(:,:,:)+(2*sigma*gfi(:,:,:,1)-sifiu*gu(:,:,:,1)).*gfi(:,:,:,1)-2*Afi*ep*gmu(:,:,:,1).*gfi(:,:,:,1) ... 
                   +2*Afi*ep*mu.*ggfi1(:,:,:,1));
    str(:,:,:,2,2)=eta*(P(:,:,:)+(2*sigma*gfi(:,:,:,2)-sifiu*gu(:,:,:,2)).*gfi(:,:,:,2)-2*Afi*ep*gmu(:,:,:,2).*gfi(:,:,:,2) ...
                   +2*Afi*ep*mu.*ggfi2(:,:,:,2));
    str(:,:,:,3,3)=eta*(P(:,:,:)+(2*sigma*gfi(:,:,:,3)-sifiu*gu(:,:,:,3)).*gfi(:,:,:,3)-2*Afi*ep*gmu(:,:,:,3).*gfi(:,:,:,3) ...
                   +2*Afi*ep*mu.*ggfi3(:,:,:,3));
 
    str(:,:,:,1,2)=eta*((2*sigma*gfi(:,:,:,2)-sifiu*gu(:,:,:,2)).*gfi(:,:,:,1)-2*Afi*ep*gmu(:,:,:,2).*gfi(:,:,:,1) ...
                   +2*Afi*ep*mu.*ggfi2(:,:,:,1));
    str(:,:,:,1,3)=eta*((2*sigma*gfi(:,:,:,3)-sifiu*gu(:,:,:,3)).*gfi(:,:,:,1)-2*Afi*ep*gmu(:,:,:,3).*gfi(:,:,:,1) ...
                   +2*Afi*ep*mu.*ggfi3(:,:,:,1));
    str(:,:,:,2,1)=eta*((2*sigma*gfi(:,:,:,1)-sifiu*gu(:,:,:,1)).*gfi(:,:,:,2)-2*Afi*ep*gmu(:,:,:,1).*gfi(:,:,:,2) ...
                   +2*Afi*ep*mu.*ggfi1(:,:,:,2));
    str(:,:,:,2,3)=eta*((2*sigma*gfi(:,:,:,3)-sifiu*gu(:,:,:,3)).*gfi(:,:,:,2)-2*Afi*ep*gmu(:,:,:,3).*gfi(:,:,:,2) ...
                   +2*Afi*ep*mu.*ggfi3(:,:,:,2));
    str(:,:,:,3,1)=eta*((2*sigma*gfi(:,:,:,1)-sifiu*gu(:,:,:,1)).*gfi(:,:,:,3)-2*Afi*ep*gmu(:,:,:,1).*gfi(:,:,:,3) ... 
                   +2*Afi*ep*mu.*ggfi1(:,:,:,3));
    str(:,:,:,3,2)=eta*((2*sigma*gfi(:,:,:,2)-sifiu*gu(:,:,:,2)).*gfi(:,:,:,3)-2*Afi*ep*gmu(:,:,:,2).*gfi(:,:,:,3) ...
                   +2*Afi*ep*mu.*ggfi2(:,:,:,3)); 
        
 
 
 dS1(:,:,:)=str(:,:,:,1,1).*gFu(:,:,:,1)+str(:,:,:,1,2).*gFu(:,:,:,2)+str(:,:,:,1,3).*gFu(:,:,:,3);
 dS2(:,:,:)=str(:,:,:,2,1).*gFu(:,:,:,1)+str(:,:,:,2,2).*gFu(:,:,:,2)+str(:,:,:,2,3).*gFu(:,:,:,3);
 dS3(:,:,:)=str(:,:,:,3,1).*gFu(:,:,:,1)+str(:,:,:,3,2).*gFu(:,:,:,2)+str(:,:,:,3,3).*gFu(:,:,:,3);
 
 gs1=grad3DR(dS1(:,:,:));
 gs2=grad3DR(dS2(:,:,:));
 gs3=grad3DR(dS3(:,:,:));
 S=gs1(:,:,:,1)+gs2(:,:,:,2)+gs3(:,:,:,3);
 
%%            dynamical equations,  for conservation of mass use  fi=fi-dt*(F+Fs);

         I=200.*u;
         %I(abs(fi)>=.9)=0;
         %I;
         %I=100*(F)*sum(sum(sum((fi>=-.99))))/Nx/Ny/Nz; % definicion de isoformasifiu.m
         %I = 120*sum(sum(sum(fi.*Gu)));
         %I = 120*Gu*sum(sum(sum((fi>=-0.99))))/Nx/Ny/Nz;
         %I(fi<=0)=0;
E=F-2*sigma*Fs+Ft+F1;
%I = 120*sum(sum(sum(E)));
	 
         

lapE=lapf3D(E);
%fi=fi-dt*(F+Fs+Ft);
fi=fi+Dfi*dt*(lapE+I);
fi(:,:,1)=fi(:,:,2);
       % Iu=sum(sum(sum(u)))/Nx/Ny/Nz-u0;
u=u+dt*(Du*lapFu+S);
        %uts
       % u(:,:,1)=u(:,:,2);  
           
    end
    ux(:,:)=u(:,Ny/2,:);
    fix(:,:)=fi(Nx/2,:,:);

    %%   
    h=isnan(fi(Nx/2,Ny/2,Nz/2));
    if h==1;
        break
    end
      

    
    %%


%para hacer peliculas 3D usar cine3D y cine2D

    Fm(:,:,:,iter+1)=fi(:,:,:);
    Um(:,:,:,iter+1)=u(:,:,:);
    %Sm(:,:,:,iter)=S(:,:,:);
   
end

time = toc(t);

save('feb3c','Fm','Um','time');




