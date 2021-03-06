function multifase(cont1,cont2,cont3)

% f = size(Fm);
Nx = 40;
Ny = 40;
Nz =70;

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

%Definiendo los valores de u y fi que serás utilizados
global Fm Um Vm
fi0 = Fm(:,:,:,1);
fi = Fm(:,:,:,cont1-1);
u = Um(:,:,:,cont1-1);
v = Vm(:,:,:,cont1-1);

disp('multifase')
for i = cont1:cont2
   for j = 1:cont3
       

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
       fi(1,:,:) = fi(2,:,:);
       fi(Nx,:,:) = fi(Nx-1,:,:);
       fi(:,1,:) = fi(:,2,:);
       fi(:,Ny,:) = fi(:,Ny-1,:);
       fi(:,:,1) = fi(:,:,2);
       fi(:,:,Nz) = fi(:,:,Nz-1);
       u(1,:,:) = u(2,:,:);
       u(Nx,:,:) = u(Nx-1,:,:);
       u(:,1,:) = u(:,2,:);
       u(:,Ny,:) = u(:,Ny-1,:);
       u(:,:,1) = u(:,:,2);
       u(:,:,Nz) = u(:,:,Nz-1);
       
       %fijando la base de la membrana
%       fi(:,:,1) = fi0(:,:,1);
%       v(1,:,:) = v(2,:,:);
%       v(Nx,:,:) = v(Nx-1,:,:);
%       v(:,1,:) = v(:,2,:);
%       v(:,Ny,:) = v(:,Ny-1,:);
%       noFlux(fi,fi);
        %noFlux2(fi,u);
        
      %condicion para parar el proceso en caso de que fi tenga entradas
      %tipo NaN

        h=isnan(fi(Nx/2,Ny/2,Nz/2));
        if h==1;
            break
        end
   end
   Fm(:,:,:,i) = fi;
   Um(:,:,:,i) = u;
   Vm(:,:,:,i) = v;
   disp(i)
end
