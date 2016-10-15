Nx = 50;
Ny = 50;
Nz = 70;


%Contadores para el proceso
cont1 = 2;
cont2 = 302;
cont3 = 200;
cont4 = 20;
cont5 = 50;

%Condicion inicial del meristemo
fi=ones(Nx,Ny,Nz);
r = zeros(Nx,Ny,Nz);

for i=1:Nx
    for j=1:Ny
        for k=1:Nz
      r(i,j,k)=sqrt((i-Nx/2)^2+(j-Ny/2)^2+(k-Nz/2)^2);
      if r(i,j,k)>=15
      fi(i,j,k)=-1;
      end
        end
   
   end
end

%Condicion inicial para el morfogeno

N = 3;
[X,Y,Z]=meshgrid(1:Nx,1:Ny,1:Nz);
teta=atan2((Y-Ny/2),(X-Nx/2));
rad=sqrt((X-Nx/2+.5).^2+(Y-Ny/2+.5).^2);
%u=1.5*exp(-((X-Nx/2).^2+(Y-Ny/2).^2+(Z-7).^2)/50);
%u=2.5*rad.*(cos(teta*N)+sin(teta*N)).*(Z/Nz)/max(max(max(rad)))+(exp(-((-X+Nx/2-.5).^2+(-Y+Ny/2-.5).^2+(-Z+R+14).^2)/80));
u=-2.5*rad.*(cos(teta*N)+sin(teta*N)).*(Z/Nz)/max(max(max(rad)))+(exp(-((X-Nx/2).^2+(Y-Ny/2).^2+(Z-1).^2)/50));
%u=2.5*rand(Nx,Ny,Nz);
%u=.1*u+.2*(rand(Nx,Ny,Nz)-.5);
%u = 1-2*u;
v=.1*u+.2*(rand(Nx,Ny,Nz)-.5);


%Definiendo Fm y Um
global Fm Um Vm
Fm = zeros(Nx,Ny,Nz,cont2+2*cont4+cont5+3);
Um = zeros(Nx,Ny,Nz,cont2+2*cont4+cont5+3);
Vm = zeros(Nx,Ny,Nz,cont2+2*cont4+cont5+3);
Fm(:,:,:,1) = fi;
Um(:,:,:,1) = u;
Vm(:,:,:,1) = v;

disp(1)

%Funcion tic toc
t = tic();

bvam2(cont1,cont2,cont3);

multifase(cont2+1,cont2+1+cont4,cont3);

bvam2(cont2+cont4+2,cont2+cont4+2+cont5,cont3);

multifase(cont2+cont4+cont5+3,cont2+2*cont4+cont5+3,cont3);

time = toc(t)/60;

save('oct14j.mat')
