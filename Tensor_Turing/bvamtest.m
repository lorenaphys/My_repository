%%% Modelo BVAM sin tensor de esfuerzos

Nx = 40;
Ny = 40;
Nz = 70;
NF = 200;
step = 50;
%dt = 1e-5;
dt = 0.0005;

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
Um = zeros(Nx,Ny,Nz,NF+1);
Vm = zeros(Nx,Ny,Nz,NF+1);

Um(:,:,:,1) = u;
Vm(:,:,:,1) = v;
M = struct('cdata',[],'colormap',[]);

%funcion que contabiliza el timepo de proceso
t = tic();

 for i = 1:NF
    for j = 1:step

      %implementacion del modelo bvam
      %u=.1*u+.2*(rand(Nx,Ny,Nz)-.5);
      
      lapu = lapf3D(u);
      lapv = lapf3D(v);
      %u = u + dt*(D*lapu +eta1*(u+a*v-C*u.*v-u.*v.^2));
      %v = v + dt*(lapv +eta1*(b*v+h*u+C*u.*v+u.*v.^2));
      
      u = u + dt*(Du*lapu + u+a*v-C*u.*v-u.*v.^2);
      v = v + dt*(Dv*lapv + b*v+h*u+C*u.*v+u.*v.^2);
      
      %condiciones de frontera
      u(:,:,1) = u(:,:,2);
      v(:,:,1) = v(:,:,2);
      
      
    end
    Um(:,:,:,i+1) = u;
    Vm(:,:,:,i+1) = v;
%     R = 11;
%     [x,y,z] = meshgrid(1:1:Ny,1:1:Nx,1:1:Nz);
%     xslice = [Nx/2-R:R:Nx/2+R,Nx/2-R:R:Nx/2+R];yslice = [Ny/2-R:R:Ny/2+R,Ny/2-R:R:Ny/2+R]; zslice = [0:3:R,0:3:R];
%     p3=slice(x,y,z,u,xslice,yslice,zslice);
%     set(p3,'FaceColor','flat','EdgeColor','none','FaceAlpha',0.1);
%     rs=max(abs(max(max(max(u)))),abs(min(min(min(u)))));
%     axis equal, view(74,18), 
%     set(gca,[-rs,rs])
%     colormap hsv;
%     N(iter) = getframe;
end
time = toc(t);


for k = 1:NF+1
    R = 11;
    u = Um(:,:,:,k);
%     [x,y,z] = meshgrid(1:1:Ny,1:1:Nx,1:1:Nz);
%     xslice = [Nx/2-R:R:Nx/2+R,Nx/2-R:R:Nx/2+R];yslice = [Ny/2-R:R:Ny/2+R,Ny/2-R:R:Ny/2+R]; zslice = [0:3:R,0:3:R];
%     p3=slice(x,y,z,u,xslice,yslice,zslice);
%     set(p3,'FaceColor','flat','EdgeColor','none','FaceAlpha',0.1);
%     rs=max(abs(max(max(max(u)))),abs(min(min(min(u)))));
%     axis equal, view(74,18), 
%     set(gca,'CLim',[-rs,rs])
%     colormap hsv;

    cdata = smooth3((Um(:,:,:,k)-min(min(min(Um(:,:,:,k)))))./...
            (max(max(max(Um(:,:,:,k))))-min(min(min(Um(:,:,:,k))))),'box',5);
    [x,y,z] = meshgrid(1:1:Nx,1:1:Ny,1:1:Nz);
    xslice = [0,Nx];yslice = [0,Ny]; zslice = [0,Nz];
    p3=slice(x,y,z,Um(:,:,:,k),xslice,yslice,zslice);
    set(p3,'FaceColor','interp','EdgeColor','none','FaceAlpha',0.5),
    axis equal, view(-70,20)
    colormap jet,
    M(k) = getframe;
end

save('junio10a');                 