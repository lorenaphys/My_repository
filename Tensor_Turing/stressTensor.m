%%%Modelo con el tensor de esfuerzos en el que se asocia al 
%%% modelo de multifases mediante su divergencia

%Definiendo las variables del modelo
Nx = 40;
Ny = 40;
NF = 20;
step = 10;
Nz = 70;
ep = 2;
Afi = 0.5;
As = 2;
Af = 2;
beta = 0.5;
sigma = 0.1;
u1 = 0;
u2 = 1;
u3 = 0;
Dfi = 1;
Du = 1;
L = 0.1;
alpha = 120;
dt = 1e-5;

%Condicion inicial para fi
fi = ones(Nx,Ny,Nz);
r = zeros(Nx,Ny,Nz);

for i = 1:Nx
	for j = 1:Ny
		for k = 1:Nz
			r(i,j,k) = sqrt((i-Nx/2)^2 + (j - Ny)^2 + (k)^2);
			if r(i,j,k) >= 15
				fi(i,j,k) = -1;
			end
		end
	end
end

%Condicion inicial para u
[x,y,z] = meshgrid(1:Nx,1:Ny,1:Nz)
%theta = atan2((y-Ny/2),(x-Nx/2));
%rad = sqrt((x-Nx/2-0.5)^2 + (y -Ny/2-0.5)^2);
u = 1.5*exp(((x-Nx/2)^2 + (y-Ny/2)^2 + (z -7)^2)/50);

%Variables que guardan a fi y u en cada iteracion
Fm = zeros(Nx,Ny,Nz,NF+1);
Um = zeros(Nx,Ny,Nz,NF+1);
Fm(:,:,:,1) = fi;
Um(:,:,:,1) = u;

%Restringiendo u a que solo actue en la membrana
u(fi<=-0.99) = 0;

%Funcion que contabiliza el tiempo del proceso
t = tic() 

%Integracion numerica
for i = 1:NF
	for j = 1:step
			lapfi = lapf3D(fi);

			%Potencial quimico
			mu = (fi.^2-1).*(fi-ep*beta*u) - ep^2*lapfi;

			%Variacion de la energia libre con respecto a fi
			lapmu = lapf3D(mu);
			varFfi = 2*Afi*mu.*(3*fi.^2-1-2*ep*beta*fi.*u) + 4*As*fi.*(fi.^2-1).*(u-u1).^2.*(u-u2).^2 + 2*Af*fi.*(u-u3)^2 - ...
			         2*sigma*lapfi - 2*Afi*ep^2*lapmu;

			%Variacion de la energia libre con respecto a u
			lapu = lapf3D(u);
			varFu = -2*Afi*ep*beta*mu.*(fi.^2-1) + 2*As*(fi.^2-1).^2.*(u-u1).*(u-u2).*(2*u-u1-u2) + 2*Af*fi.^2.*(u-u3) - ...
			        As*L*lapu;
	end
end
