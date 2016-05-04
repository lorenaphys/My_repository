%%%Modelo con el tensor de esfuerzos en el que se asocia al 
%%% modelo de multifases mediante su divergencia

%Definiendo las variables del modelo
Nx = 40;
Ny = 40;
Nz = 70;
NF = 100;
step = 30;
ep = 1e-5;
Afi = 0.5;
As = 0.05;
Af = 0.05;
beta = 0.1;
sigma = -0.1;
u1 = 0;
u2 = 1;
u3 = 0;
Dfi = 1;
Du = 50;
L = -0.25;
alpha = 0.01;
dt = 1e-5;
%eta = 1;

%Condicion inicial para fi
fi = ones(Nx,Ny,Nz);
r = zeros(Nx,Ny,Nz);

for i = 1:Nx
	for j = 1:Ny
		for k = 1:Nz
			r(i,j,k) = sqrt((i-Nx/2)^2 + (j - Ny/2)^2 + (k)^2);
			if r(i,j,k) >= 15
				fi(i,j,k) = -1;
			end
		end
	end
end

%Condicion inicial para u
[x,y,z] = meshgrid(1:Nx,1:Ny,1:Nz);
%theta = atan2((y-Ny/2),(x-Nx/2));
%rad = sqrt((x-Nx/2-0.5)^2 + (y -Ny/2-0.5)^2);
u=1.5*exp(-((x-Nx/3-.5).^2+(y-Ny/3-.5).^2+(z-13.5).^2)/50);
%u=2.5*rad.*(cos(teta*N)+sin(teta*N)).*(Z/Nz)/max(max(max(rad)))+(exp(-((-X+Nx/2-.5).^2+(-Y+Ny/2-.5).^2+(-Z+R+14).^2)/80));
%u=-2.5*rad.*(cos(teta*N)+sin(teta*N)).*(Z/Nz)/max(max(max(rad)))+(exp(-((X-Nx/2).^2+(Y-Ny/2).^2+(Z-R+14).^2)/50));
%u=2.5*rand(Nx,Ny,Nz);

%Variables que guardan a fi y u en cada iteracion
Fm = zeros(Nx,Ny,Nz,NF+1);
Um = zeros(Nx,Ny,Nz,NF+1);
Fm(:,:,:,1) = fi;
Um(:,:,:,1) = u;

%Restringiendo u a que solo actue en la membrana
u(fi<=-0.99) = 0;

%Funcion que contabiliza el tiempo del proceso
t = tic(); 

%Integracion numerica
for i = 1:NF 
    for j = 1:step
			lapfi = lapf3D(fi);

			%Potencial quimico
			mu = (fi.^2-1).*(fi-ep*beta*u) - ep^2*lapfi;

			%Variacion de la energia libre con respecto a fi
			lapmu = lapf3D(mu);
			varFfi = 2*Afi*mu.*(3*fi.^2-1-2*ep*beta*fi.*u) + 4*As*fi.*(fi.^2-1).*(u-u1).^2.*(u-u2).^2 + 2*Af*fi.*(u-u3).^2 - ...
			         2*sigma*lapfi - 2*Afi*ep^2*lapmu;

			%Variacion de la energia libre con respecto a u
			lapu = lapf3D(u);
			varFu = -2*Afi*ep*beta*mu.*(fi.^2-1) + 2*As*(fi.^2-1).^2.*(u-u1).*(u-u2).*(2*u-u1-u2) + 2*Af*fi.^2.*(u-u3) - ...
			        As*L*lapu;
                
            %Tensor de esfuerzos
            gfi = grad3DR(fi);
            gmu = grad3DR(mu);
            gu = grad3DR(u);
            Vs = (fi.^2-1).^2.*(u-u1).^2.*(u-u2).^2 + L*abs(gu(:,:,:,1).^2+gu(:,:,:,2).^2+gu(:,:,:,3).^2);
            Vf = fi.^2.*(u-u3).^2;
            P = Afi*mu.^2 + As*Vs + Af*Vf + sigma*abs(gfi(:,:,:,1).^2+gfi(:,:,:,2).^2+gfi(:,:,:,3).^2) - fi.*varFfi;
            %P=(Afi*mu.^2-sigma*abs(gfi(:,:,:,1).^2+gfi(:,:,:,2).^2+gfi(:,:,:,3).^2)+abs(gu(:,:,:,1).^2+gu(:,:,:,2).^2+gu(:,:,:,3).^2)...
%               +sifiu*(gfi(:,:,:,1).*gu(:,:,:,1)+gfi(:,:,:,2).*gu(:,:,:,2)+gfi(:,:,:,3).*gu(:,:,:,3)) +As*(fi.^2-1).^2.*(u-u1).^2.*(u-u2).^2 ...
%               +Af*fi.^2.*(u-u3).^2-fi(:,:,:).*(F(:,:,:)-2*sigma*Fs(:,:,:)+Ft(:,:,:)+F1(:,:,:)));
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
                         
            %Definiendo la relacion entre varFfi y el gradiente del tensor
            gst11 = grad3DR(str(:,:,:,1,1));
            gst22 = grad3DR(str(:,:,:,2,2));
            gst33 = grad3DR(str(:,:,:,3,3));
            
            gst12 = grad3DR(str(:,:,:,1,2));
            gst13 = grad3DR(str(:,:,:,1,3));
            gst21 = grad3DR(str(:,:,:,2,1));
            gst23 = grad3DR(str(:,:,:,2,3));
            gst31 = grad3DR(str(:,:,:,3,1));
            gst32 = grad3DR(str(:,:,:,3,2));
            
            ggst11 = grad3DR(gst11(:,:,:,1)./fi);
            ggst22 = grad3DR(gst11(:,:,:,2)./fi);
            ggst33 = grad3DR(gst11(:,:,:,3)./fi);
            
            ggst12 = grad3DR(gst11(:,:,:,2)./fi);
            ggst13 = grad3DR(gst11(:,:,:,3)./fi);
            ggst21 = grad3DR(gst11(:,:,:,1)./fi);
            ggst23 = grad3DR(gst11(:,:,:,3)./fi);
            ggst31 = grad3DR(gst11(:,:,:,1)./fi);
            ggst32 = grad3DR(gst11(:,:,:,2)./fi);
            lapFfi = -(ggst11(:,:,:,1) + ggst22(:,:,:,2) + ggst33(:,:,:,3) + ggst12(:,:,:,1) + ggst13(:,:,:,1) + ggst21(:,:,:,2) +...
                     ggst23(:,:,:,2) + ggst31(:,:,:,3) + ggst32(:,:,:,3));
               
            %Crecimiento de fi debido a la sustancia
            lapFu = lapf3D(varFu);
            I = lapFu*sum(sum(sum(fi >= -0.99)));
      
            %Dinamica de los campos fase
            fi = fi + Dfi*dt*(lapFfi + alpha*I);
            u = u + Du*dt*lapFu;
            
            %condiciones de frontera
%             fi(1,:,:) = fi(2,:,:);
%             fi(Nx,:,:) = fi(Nx-1,:,:);
%             fi(:,1,:) = fi(:,2,:);
%             fi(:,Ny,:) = fi(:,Ny-1,:);
%             u(1,:,:) = u(2,:,:);
%             u(Nx,:,:) = u(Nx-1,:,:);
%             u(:,1,:) = u(:,2,:);
%             u(:,Ny,:) = u(:,Ny-1,:);

              noFlux(fi,fi);
              noFlux(fi,u);

            %condicion para parar el proceso en caso de que fi tenga entradas
            %tipo NaN
            k=isnan(fi(Nx/2,Ny/2,Nz/2));
            if k == 1;
                break
            end
    end
   Fm(:,:,:,i+1) = fi;
   Um(:,:,:,i+1) = u;
end

time = toc(t);
save('mayo3d')
