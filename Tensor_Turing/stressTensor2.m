%% Segunda opcion para la energia libre

Nx = 40;
Ny = 40;
Nz = 70;
NF = 100;
step = 5;
ep = 2;
Afi = 0.5;
As = 2;
Af = 2;
Dfi = 1;
Du = 10;
sigma = -0.1;
beta = 0.1;
L = 0.1;
eta = 5;
alpha = 200;
dt = 1e-5;

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

%R = 9;
%N = 6;
[X,Y,Z]=meshgrid(1:Nx,1:Ny,1:Nz);
%teta=atan2((Y-Ny/2),(X-Nx/2));
%rad=sqrt((X-Nx/2+.5).^2+(Y-Ny/2+.5).^2);
u=1.5*exp(-((X-Nx/2-.5).^2+(Y-Ny/2-.5).^2+(Z-7).^2)/50);
%u=2.5*rad.*(cos(teta*N)+sin(teta*N)).*(Z/Nz)/max(max(max(rad)))+(exp(-((-X+Nx/2-.5).^2+(-Y+Ny/2-.5).^2+(-Z+R+14).^2)/80));
%u=-2.5*rad.*(cos(teta*N)+sin(teta*N)).*(Z/Nz)/max(max(max(rad)))+(exp(-((X-Nx/2).^2+(Y-Ny/2).^2+(Z-R+14).^2)/50));
%u=2.5*rand(Nx,Ny,Nz);

Fm = zeros(Nx,Ny,Nz,NF+1);
Um = zeros(Nx,Ny,Nz,NF+1);

Fm(:,:,:,1) = fi;
Um(:,:,:,1) = u;

u(fi<=-0.99) = 0;

t = tic();
for i = 1:NF
   for j = 1:step
      lapfi = lapf3D(fi);
      mu = (fi.^2-1).*(fi-ep*beta*u)-ep.^2*lapfi;
      lapmu = lapf3D(mu);
      
      varFfi = 2*Afi*mu.*(3*fi.^2-1-2*ep*beta*fi.*u)-2*sigma*lapfi...
               -2*Afi*ep^2*lapmu;
      
      lapu = lapf3D(u);
      varFu = -2*Afi*ep*beta*mu.*(fi.^2-1)-2*L*lapu;
      
      %%   tensor de esfuerzos
      gfi=grad3DR(fi);
      gmu=grad3DR(mu);
      gu=grad3DR(u);
            
      P = Afi*mu.^2 + L*(gu(:,:,:,1).^2+gu(:,:,:,2).^2+gu(:,:,:,3).^2) +...
          sigma*abs(gfi(:,:,:,1).^2+gfi(:,:,:,2).^2+gfi(:,:,:,3).^2) - ...
          fi.*varFfi;
           
      ggfi1 = grad3DR(gfi(:,:,:,1));
      ggfi2 = grad3DR(gfi(:,:,:,2));
      ggfi3 = grad3DR(gfi(:,:,:,3));        
        

      str=zeros(Nx,Ny,Nz,3,3);
      str(:,:,:,1,1) = P(:,:,:)-2*sigma*gfi(:,:,:,1).*gfi(:,:,:,1)-2*Afi*ep*gmu(:,:,:,1).*gfi(:,:,:,1) ... 
                       +2*Afi*ep*mu.*ggfi1(:,:,:,1);
      str(:,:,:,2,2) = P(:,:,:)-2*sigma*gfi(:,:,:,2).*gfi(:,:,:,2)-2*Afi*ep*gmu(:,:,:,2).*gfi(:,:,:,2) ...
                       +2*Afi*ep*mu.*ggfi2(:,:,:,2);
      str(:,:,:,3,3) = P(:,:,:)-2*sigma*gfi(:,:,:,3).*gfi(:,:,:,3)-2*Afi*ep*gmu(:,:,:,3).*gfi(:,:,:,3) ...
                       +2*Afi*ep*mu.*ggfi3(:,:,:,3);
 
      str(:,:,:,1,2) = -2*sigma*gfi(:,:,:,2).*gfi(:,:,:,1)-2*Afi*ep*gmu(:,:,:,2).*gfi(:,:,:,1) ...
                       +2*Afi*ep*mu.*ggfi2(:,:,:,1);
      str(:,:,:,1,3) = -2*sigma*gfi(:,:,:,3).*gfi(:,:,:,1)-2*Afi*ep*gmu(:,:,:,3).*gfi(:,:,:,1) ...
                       +2*Afi*ep*mu.*ggfi3(:,:,:,1);
      str(:,:,:,2,1) = -2*sigma*gfi(:,:,:,1).*gfi(:,:,:,2)-2*Afi*ep*gmu(:,:,:,1).*gfi(:,:,:,2) ...
                       +2*Afi*ep*mu.*ggfi1(:,:,:,2);
      str(:,:,:,2,3) = -2*sigma*gfi(:,:,:,3).*gfi(:,:,:,2)-2*Afi*ep*gmu(:,:,:,3).*gfi(:,:,:,2) ...
                       +2*Afi*ep*mu.*ggfi3(:,:,:,2);
      str(:,:,:,3,1) = -2*sigma*gfi(:,:,:,1).*gfi(:,:,:,3)-2*Afi*ep*gmu(:,:,:,1).*gfi(:,:,:,3) ... 
                       +2*Afi*ep*mu.*ggfi1(:,:,:,3);
      str(:,:,:,3,2) = -2*sigma*gfi(:,:,:,2).*gfi(:,:,:,3)-2*Afi*ep*gmu(:,:,:,2).*gfi(:,:,:,3) ...
                       +2*Afi*ep*mu.*ggfi2(:,:,:,3); 
        
 
      gFu=grad3DR(varFu);
      dS1(:,:,:)=str(:,:,:,1,1).*gFu(:,:,:,1)+str(:,:,:,1,2).*gFu(:,:,:,2)+str(:,:,:,1,3).*gFu(:,:,:,3);
      dS2(:,:,:)=str(:,:,:,2,1).*gFu(:,:,:,1)+str(:,:,:,2,2).*gFu(:,:,:,2)+str(:,:,:,2,3).*gFu(:,:,:,3);
      dS3(:,:,:)=str(:,:,:,3,1).*gFu(:,:,:,1)+str(:,:,:,3,2).*gFu(:,:,:,2)+str(:,:,:,3,3).*gFu(:,:,:,3);
 
      gs1=grad3DR(dS1(:,:,:));
      gs2=grad3DR(dS2(:,:,:));
      gs3=grad3DR(dS3(:,:,:));
      S=gs1(:,:,:,1)+gs2(:,:,:,2)+gs3(:,:,:,3);
      
      %%ecuaciones dinamicas
      %m = varFu.*sum(sum(sum(fi>=-0.99)));
      m = u;
      
      lapFfi = lapf3D(varFfi);
      fi = fi+Dfi*dt*(lapFfi+alpha*m);
      
      lapFu = lapf3D(varFu);
      u = u+Du*dt*(lapFu+eta*S);
      %%condiciones de frontera
      noFlux2(fi,u);
      
      h=isnan(fi(Nx/2,Ny/2,Nz/2));
      if h==1;
          break
      end
   end
   Fm(:,:,:,NF+1) = fi;
   Um(:,:,:,NF+1) = u;
end

time = toc(t);
save('junio21c');