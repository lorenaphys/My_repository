%%% Prueba de cuańto valen los términos involucrados en las dinámicas de fi
%%% y u para así determinar qué valor deben tener las constantes que
%%% indican la importancia que tienen en el modelo
ep1=2;
ep=ep1^2;
sigma=-0.1;
Nx=40;
Ny=40;
Nz=70;
sifiu=0.; 
Du=10;
Dfi=1;
eta=0.1;
Afi = 0.5;
Af = 1;
As = 2; 
beta = 0.1;
u1 = 0;
u2 = 1;
u3 = 0;
L = 0.1;

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

u = zeros(Nx,Ny,Nz);
        
     for i=1:Nx
         for j=1:Ny
             for k=1:Nz
                 u(i,j,k)=1.5*exp(-((i-Nx/2-0.5)^2+(j-Ny/2-0.5)^2+(k-7-0.5)^2)/50);%primera prueba
             end
         end
     end
     
dt=1e-5;
     
lapfi = lapf3D(fi);
        
mu=(fi-ep1*beta*u).*((fi).^2-1)-ep*lapfi;        
        
lapmu = lapf3D(mu);
        
F=Afi*((3*fi.^2-1-2*ep1*beta*fi.*u).*mu-2*ep*lapmu);     
        
lapmu = lapf3D(mu);

lapu = lapf3D(u);
        
Fs=lapfi;
        
%lapFs = lapf3D(Fs);

Gup = 2*As*(fi.^2-1).*(u-u1).*(u-u2).*(2*u-u1-u2)+2*Af*fi.^2.*(u-u3);
        Gu=-2*Afi*ep1*beta*mu.*((fi).^2-1)-sifiu*lapfi+Gup; 

       
        Fu=-2*As*L*lapu+Gup+Gu;
        gFu=grad3DR(Fu);
       
	lapFu = lapf3D(Fu);
        
        Ft=-sifiu*lapu; %  surface tension between membrane and u
        
	%lapFt = lapf3D(Ft);
        
        F1 = 4*As*fi.*(fi.^2-1).*(u-u1).^2.*(u-u2).^2 +2*Af*fi.*(u-u3).^2;
        
	%lapF1 = lapf3D(F1);
        
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
 
 I=200.*u;
 
 %I(abs(fi)>=0.9)=0;
E=F-2*sigma*Fs+Ft+F1;

lapE=lapf3D(E);

fi=fi+Dfi*dt*(lapE+I);
fi(:,:,1)=fi(:,:,2);

u=u+dt*(Du*lapFu+S);

u(:,:,1)=u(:,:,2);
 