function bvam2(U,Fi,cont1,cont2,cont3)
%%% cont1 indica con qué número inicia el primer loop
%%% cont2 indica con qué número termina el primer loop
%%% cont3 indica con qué número termina el segundo loop, ya que este siempre inicia en 1

f = size(Fi
U = u;
Nx = f(1);
Ny = f(2);
Nz = f(3);

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

disp('bvam')
for i = cont1:cont2
	for j = 1:cont3
        	lapu = lapf3D(u);
        	lapv = lapf3D(v);
        	u = u + dt1*(Du*lapu + u+a*v-C*u.*v-u.*v.^2);
        	v = v + dt1*(Dv*lapv + b*v+h*u+C*u.*v+u.*v.^2);
	end
	Um(:,:,:,i) = u;
	Fm(:,:,:,i) = Fi;
	disp(i)
end
