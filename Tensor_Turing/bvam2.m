function bvam2(cont1,cont2,cont3)
%%% cont1 indica con qué número inicia el primer loop
%%% cont2 indica con qué número termina el primer loop
%%% cont3 indica con qué número termina el segundo loop, ya que este siempre inicia en 1

% f = size(Fm);
 %Nx = 40;
 %Ny = 40;
 %Nz = 70;
dt = 1e-5;
dt1 = 500*dt;

%Parametros del modelo BVAM

h = -1;
C = 1.57;

%Primer conjunto, para kc = 0.46 (ac = 1.121)

eta = sqrt(3);
%eta = 3;
D = 0.516;
Du = D/eta;
Dv = 1/eta;
a = 1/0.899;
b = -0.91/0.899;
% eta1 = 0.450;

%Segundo conujunto, para kc = 0.85 (ac = 2.583)

% D = 0.122;
% a = 2.513;
% b = -1.005;
% eta1 = 0.199;

%Definiendo los valores de u, y v
global Fm Um Vm
u = Um(:,:,:,cont1-1);
v = Vm(:,:,:,cont1-1);
fi = Fm(:,:,:,cont1-1);
f = size(Fm);
Nx = f(1);
Ny = f(2);
Nz = f(3);

disp('BVAM')
for i = cont1:cont2
	for j = 1:cont3
        	lapu = lapf3D(u);
        	lapv = lapf3D(v);
        	u = u + dt1*(Du*lapu + u+a*v-C*u.*v-u.*v.^2);
        	v = v + dt1*(Dv*lapv + b*v+h*u+C*u.*v+u.*v.^2);
        	%condiciones de frontera
			u(1,:,:) = u(2,:,:);
			u(Nx,:,:) = u(Nx-1,:,:);
			u(:,1,:) = u(:,2,:);
			u(:,Ny,:) = u(:,Ny-1,:);
			u(:,:,1) = u(:,:,2);
			u(:,:,Nz) = u(:,:,Nz-1);
			v(1,:,:) = v(2,:,:);
			v(Nx,:,:) = v(Nx-1,:,:);
			v(:,1,:) = v(:,2,:);
			v(:,Ny,:) = v(:,Ny-1,:);
			v(:,:,1) = v(:,:,2);
			v(:,:,Nz) = v(:,:,Nz-1);
            %u(fi<=-0.99) = 0;
            %v(fi<=-0.99) = 0;
	end
	Um(:,:,:,i) = u;
	Vm(:,:,:,i) = v;
	Fm(:,:,:,i) = fi;
	u(fi<=-0.99) = 0;
	v(fi<=-0.99) = 0;
	disp(i)
end


