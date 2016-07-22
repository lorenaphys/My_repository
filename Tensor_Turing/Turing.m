Nx = 40;
Ny = 40;
Nz = 70;
Afi = 1;
As = 0.01;
Af = 0.01;
sigma = -0.1;
ep = 1;
Dfi = 1;
du =1;
u1 = 0;
u2 = 1;
u3 = 0;
beta = 0.5;
L = 0.07;
alpha = 120;
dt = 1e-5;
dt1 = 500*dt;

%Condicion inicial del meristemo
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

%Condicion inicial para el morfogeno

N = 0;
[X,Y,Z]=meshgrid(1:Nx,1:Ny,1:Nz);
teta=atan2((Y-Ny/2),(X-Nx/2));
rad=sqrt((X-Nx/2+.5).^2+(Y-Ny/2+.5).^2);
u=1.5*exp(-((X-Nx/2-.5).^2+(Y-Ny/2-.5).^2+(Z-7).^2)/20);
%u=2.5*rad.*(cos(teta*N)+sin(teta*N)).*(Z/Nz)/max(max(max(rad)))+(exp(-((-X+Nx/2-.5).^2+(-Y+Ny/2-.5).^2+(-Z+R+14).^2)/80));
%u=-2.5*rad.*(cos(teta*N)+sin(teta*N)).*(Z/Nz)/max(max(max(rad)))+(exp(-((X-Nx/2).^2+(Y-Ny/2).^2+(Z-1).^2)/50));
%u=2.5*rand(Nx,Ny,Nz);
%u=.1*u+.2*(rand(Nx,Ny,Nz)-.5);
v=.1*u+.2*(rand(Nx,Ny,Nz)-.5);


%Definiendo Fm y Um
Fm = zeros(Nx,Ny,Nz,cont2+2*cont4+cont5+1);
Um = zeros(Nx,Ny,Nz,cont2+2*cont4+cont5+1);
Fm(:,:,:,1) = fi;
Um(:,:,:,1) = u;

disp(1)

%Contadores para el proceso
cont1 = 2;
cont2 = 351;
cont3 = 200;
cont4 = 50;
cont5 = 20;

bvam2(u,fi,cont1,cont2,cont3);

Um(:,:,:,cont2) = Ucont2;
Fm(:,:,:,cont2) = Fcont2;
multifase(Ucont2,Fi,cont1,cont2,cont3)
