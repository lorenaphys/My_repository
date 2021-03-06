%%% Modelo BVAM sin tensor de esfuerzos

Nx = 40;
Ny = 40;
Nz = 70;
NF = 20;
step = 200;
Afi = 1;
As = 0.01;%0.25
Af = 0.01;
sigma = -0.1;
ep = 1;%0.06
%Du = 1e-4;
Dfi = 1;%1
du =1;
%eta = 1;
u1 = 0;
u2 = 1;
u3 = 0;
beta = 0.5;
L = 0.07;
alpha = 120;%5
dt = 1e-5;
dt1 = 500*dt;
cont = 350;
cont2 = 200;
cont3 = 80;
cont4 = 200;

%Parametros del modelo BVAM

h = -1;
C = 0;

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

%filotaxia
%tetar=0;   % X rotation
%fir=0;     % Z rotation
%RX=(X)*cos(fir)-(Y)*sin(fir)*cos(tetar)+(Z)*sin(fir)*sin(tetar);
%RY=(X)*sin(fir)+(Y)*cos(fir)*cos(tetar)-(Z)*cos(fir)*sin(tetar);
%RZ=(Y)*sin(tetar)+(Z)*cos(tetar);
%teta=atan2((RY-Ny/2),(RX-Nx/2));
%rad=sqrt((RX-Nx/2+.25).^2+(RY-Ny/2+.25).^2);
%u=2.5*rad.*(cos(teta*N)+sin(teta*N)).*(RZ/Nz)/max(max(max(rad)))+...
%  (exp(-((RX-Nx/2).^2+(RY-Ny/2).^2+(RZ-7).^2)/50));
v=.1*u+.2*(rand(Nx,Ny,Nz)-.5);

%iteraciones del modelo
    %Variables que guardan todas la iteraciones
Fm = zeros(Nx,Ny,Nz,NF+1);
Um = zeros(Nx,Ny,Nz,NF+1+cont+cont3);
Um(:,:,:,1) = u;

%funcion que contabiliza el tiempo de proceso
t = tic();
disp('cont =')

for i = 1:cont
	for j = 1:cont2
        	lapu = lapf3D(u);
        	lapv = lapf3D(v);
        	u = u + dt1*(Du*lapu + u+a*v-C*u.*v-u.*v.^2);
        	v = v + dt1*(Dv*lapv + b*v+h*u+C*u.*v+u.*v.^2);
	end
	Um(:,:,:,i+1) = u;
	disp(i)
end

Fm(:,:,:,1) = fi;
u(fi<=-0.99) = 0;
v(fi<=-0.99) = 0;

disp('NF =')
for i = 1:NF
   for j = 1:step
       

      lapfi = lapf3D(fi); 
      %potencial quimico
      mu = (fi.^2-1).*(fi-ep*beta*u)-ep^2*lapfi;
      
      %variacion de la energia libre con respecto a fi
      lapmu = lapf3D(mu);
      varFfi = 2*Afi*mu.*(3*fi.^2-1-2*ep*beta*fi.*u) + 4*As*fi.*(fi.^2-1).*(u-u1).^2.*(u-u2).^2 + 2*Af*fi.*(u-u3).^2 - ...
               2*sigma*lapfi - 2*Afi*ep^2*lapmu;
           
      %variacion de la energia libre con respecto a u
      lapu = lapf3D(u);
      varFu = -2*Afi*ep*beta*mu.*(fi.^2-1) + 2*As*(fi.^2-1).^2.*(u-u1).*(u-u2).*(2*u-u1-u2) + 2*Af*fi.^2.*(u-u3) -...
              As*L*lapu;
           
      %crecimiento de fi debido a la sustancia
      lapFu = lapf3D(varFu);
%      I = lapFu*sum(sum(sum(fi >= -0.99)));
      I = u;      
     
      %dinamica del meristemo 
      lapFfi = lapf3D(varFfi);
      fi = fi + Dfi*dt*(lapFfi + alpha*I);
      
      %dinamica del morfogeno
      u = u + du*dt*lapFu;
      
%       u(fi<=-0.99) = 0;
%       v(fi<=-0.99) = 0;
      
      %condiciones de frontera
%       fi(1,:,:) = fi(2,:,:);
%       fi(Nx,:,:) = fi(Nx-1,:,:);
%       fi(:,1,:) = fi(:,2,:);
%       fi(:,Ny,:) = fi(:,Ny-1,:);
%       u(1,:,:) = u(2,:,:);
%       u(Nx,:,:) = u(Nx-1,:,:);
%       u(:,1,:) = u(:,2,:);
%       u(:,Ny,:) = u(:,Ny-1,:);
%       v(1,:,:) = v(2,:,:);
%       v(Nx,:,:) = v(Nx-1,:,:);
%       v(:,1,:) = v(:,2,:);
%       v(:,Ny,:) = v(:,Ny-1,:);
%       noFlux(fi,fi);
        noFlux2(fi,u);
        noFlux2(fi,v);
        
      %condicion para parar el proceso en caso de que fi tenga entradas
      %tipo NaN

        h=isnan(fi(Nx/2,Ny/2,Nz/2));
        if h==1;
            break
        end
   end
   Fm(:,:,:,i+1) = fi;
   Um(:,:,:,cont+1+i) = u;
   %implementacion del modelo bvam
   %u(fi<=-0.99) = 0;
   %v(fi<=-0.99) = 0;
   disp(i)
end

disp('cont3 = ')
%u=2.5*(rand(Nx,Ny,Nz)-.5);
%v=.2*u+0.2*(rand(Nx,Ny,Nz)-.5);
for k = 1:cont3
	for l = 1:cont4
		lapu = lapf3D(u);
		lapv = lapf3D(v);
		u = u + dt1*(Du*lapu + u+a*v-C*u.*v-u.*v.^2);
		v = v + dt1*(Dv*lapv + b*v+h*u+C*u.*v+u.*v.^2);
        	h=isnan(u(Nx/2,Ny/2,Nz/2));
        	if h==1;
            		break
        	end
	end
Um(:,:,:,cont+NF+1+k) = u;
disp(k)
end
time = toc(t);

save('julio21a');                                                                                                                               
