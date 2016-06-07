%%% Modelo BVAM sin tensor de esfuerzos

Nx = 40;
Ny = 40;
Nz = 70;
% NF = 20;
% step = 10;
dt = 1e-5;
%dt1 = 0.02;

%Parametros del modelo BVAM

h = -1;
C = 0.02;

%Primer conjunto, para kc = 0.46 (ac = 1.121)

eta = sqrt(3);
D = 0.516;
Du = D/eta;
Dv = 1/eta;
a = 1/0.899;
b = -0.91/0.899;
eta1 = 0.450;

%Segundo conujunto, para kc = 0.85 (ac = 2.583)

% D = 0.122;
% a = 2.513;
% b = -1.005;
% eta1 = 0.199;

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

 for i = 1:NF
    for j = 1:step

      %implementacion del modelo bvam
      %u=.1*u+.2*(rand(Nx,Ny,Nz)-.5);
      
      lapu = lapf3D(u);
      lapv = lapf3D(v);
      u = u + dt*(D*lapu +eta1*(u+a*v-C*u.*v-u.*v.^2));
      v = v + dt*(lapv +eta1*(b*v+h*u+C*u.*v+u.*v.^2));
      
      %u = u + dt*(Du*lapu + u+a*v-C*u.*v-u.*v.^2);
      %v = v + dt*(Dv*lapv + b*v+h*u+C*u.*v+u.*v.^2);
      
      %condiciones de frontera
      fi(:,:,1) = fi(:,:,2);
      u(:,:,1) = u(:,:,2);
      v(:,:,1) = v(:,:,2);
      
      
   end
   
end
time = toc(t);

                 