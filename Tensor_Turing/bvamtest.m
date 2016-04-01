%%% Modelo BVAM sin tensor de esfuerzos

Nx = 40;
Ny = 40;
Nz = 70;
% NF = 20;
% step = 10;
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
alpha = 1e-2;
dt = 1e-5;
%dt1 = 0.02;

%Parametros del modelo BVAM

h = -1;
C = 0;

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
      r(i,j,k)=sqrt((i-Nx/2)^2+(j-Ny/2)^2+(k)^2);
      if r(i,j,k)>=15
      fi(i,j,k)=-1;
      end
        end
   
   end
end

%Condicion inicial para el morfogeno

%N = 3;
[X,Y,Z]=meshgrid(1:Nx,1:Ny,1:Nz);
%teta=atan2((Y-Ny/2),(X-Nx/2));
%rad=sqrt((X-Nx/2+.5).^2+(Y-Ny/2+.5).^2);
%u=1.5*exp(-((X-Nx/2-.5).^2+(Y-Ny/2-.5).^2+(Z-7).^2)/50);
%u=2.5*rad.*(cos(teta*N)+sin(teta*N)).*(Z/Nz)/max(max(max(rad)))+(exp(-((-X+Nx/2-.5).^2+(-Y+Ny/2-.5).^2+(-Z+R+14).^2)/80));
%u=-2.5*rad.*(cos(teta*N)+sin(teta*N)).*(Z/Nz)/max(max(max(rad)))+(exp(-((X-Nx/2).^2+(Y-Ny/2).^2+(Z-1).^2)/50));
u=2.5*rand(Nx,Ny,Nz);
v=.1*u+.2*(rand(Nx,Ny,Nz)-.5);

%iteraciones del modelo
    %Variables que guardan todas la iteraciones

    
    
%funcion que contabiliza el timepo de proceso
t = tic();

fi0 = fi;
u0 = u;

% for i = 1:NF
%    for j = 1:step
       

      lapfi = lapf3D(fi); 
      %potencial quimico
      mu = (fi.^2-1).*(fi-ep*beta*u)-ep^2*lapfi;
      
      %variacion de la energia libre con respecto a fi
      lapmu = lapf3D(mu);
      varFfi = 2*Afi*mu.*(3*fi.^2-1-2*ep*beta*fi.*u) + 4*As*fi.*(fi.^2-1).*(u-u1).^2.*(u-u2).^2 + 2*Af*fi.*(u-u3).^2 - ...
               2*sigma*lapfi - 2*Afi*ep^2*lapmu;
           
      %variacion de la energia libre con respecto a u
      lapu = lapf3D(u);
      varFu = -2*Afi*ep*beta*mu.*(fi.^2-1) + 2*As*(fi.^2-1).^2.*(2*u-u1-u2) + 2*Af*fi.^2.*(u-u3) - As*L*lapu;
           
      %crecimiento de fi debido a la sustancia
      lapFu = lapf3D(varFu);
      I = lapFu*sum(sum(sum(fi >= -0.99)));
      
      %dinamica del meristemo
      lapFfi = lapf3D(varFfi);
      fi = fi + Dfi*dt*(lapFfi + alpha*I);
      
      %implementacion del modelo bvam
      %u=.1*u+.2*(rand(Nx,Ny,Nz)-.5);
      
      lapu = lapf3D(u);
      lapv = lapf3D(v);
      u = u + dt*(D*lapu + eta1*(u+a*v-C*u.*v-u.*v.^2));
      v = v + dt*(lapv + eta1*(b*v+h*u+C*u.*v+u.*v.^2));
      
      %condiciones de frontera
      fi(:,:,1) = fi(:,:,2);
      u(:,:,1) = u(:,:,2);
      v(:,:,1) = v(:,:,2);
      
      
   %end
   
%end
time = toc(t);

save('marzo31a');
                 