%%% Modelo BVAM sin tensor de esfuerzos

Nx = 40;
Ny = 40;
Nz = 70;
NF = 20;
step = 10;
Afi = 0.5;
As = 0.05;
Af = 0.05;
sigma = -0.1;
ep = 0.05;
Du = 1;
Dfi = 1;
eta = 1;
u1 = 0;
u2 = 1;
u3 = 0;
beta = 0.1;
L = -10;

%Parametros del modelo BVAM

h = -1;

%Primer conjunto, para kc = 0.46 (ac = 1.121)

D = 0.516;
a = 1.112;
b = -1.01;
eta1 = 0.450;

%Segundo conujunto, para kc = 0.85 (ac = 2.583)

% D = 0.122;
% a = 2.513;
% b = -1.005;
% eta1 = 0.199;

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

%Condicion inicila para el morfogeno

N = 3;
[X,Y,Z]=meshgrid(1:Nx,1:Ny,1:Nz);
teta=atan2((Y-Ny/2),(X-Nx/2));
rad=sqrt((X-Nx/2+.5).^2+(Y-Ny/2+.5).^2);
%u=1.5*exp(-((X-Nx/2-.5).^2+(Y-Ny/2-.5).^2+(Z-R+2).^2)/100);
%u=2.5*rad.*(cos(teta*N)+sin(teta*N)).*(Z/Nz)/max(max(max(rad)))+(exp(-((-X+Nx/2-.5).^2+(-Y+Ny/2-.5).^2+(-Z+R+14).^2)/80));
u=-2.5*rad.*(cos(teta*N)+sin(teta*N)).*(Z/Nz)/max(max(max(rad)))+(exp(-((X-Nx/2).^2+(Y-Ny/2).^2+(Z-1).^2)/50));
%u=2.5*rand(Nx,Ny,Nz);
                                                                                                                               