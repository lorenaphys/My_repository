Nx = 40;
Ny = 40;
Nz = 70;
fi = ones(Nx,Ny,Nz);
r = zeros(Nx,Ny,Nz);
%u = zeros(Nx,Ny,Nz);

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
    
    %Rompimiento de simetria
    R = 9;
    %N = 6;
    [X,Y,Z]=meshgrid(1:Nx,1:Ny,1:Nz);
        %teta=atan2((Y-Ny/2),(X-Nx/2));
        %rad=sqrt((X-Nx/2+.5).^2+(Y-Ny/2+.5).^2);
        u=1.5*exp(-((X-Nx/3-.5).^2+(Y-Ny/3-5-.5).^2+(Z-R+2).^2)/50);
        %u=2.5*rad.*(cos(teta*N)+sin(teta*N)).*(Z/Nz)/max(max(max(rad)))+(exp(-((-X+Nx/2-.5).^2+(-Y+Ny/2-.5).^2+(-Z+R+14).^2)/80));
        %u=-2.5*rad.*(cos(teta*N)+sin(teta*N)).*(Z/Nz)/max(max(max(rad)))+(exp(-((X-Nx/2).^2+(Y-Ny/2).^2+(Z-R+14).^2)/50));
        %u=2.5*rand(Nx,Ny,Nz);
    
%     [~,R1]=min(abs(fi(Nx/2,Ny/2,:)));
% 
%     for i=1:Nx
%         for j=1:Ny
%             for k=1:Nz
%                 u(i,j,k)=1+3.5*exp(-((i-Nx/2)^2+(j-Ny/2)^2+(k-R1+8)^2)/200);%primera prueba
%             end
%         end
%     end
    
    u(fi<=-.99)=0;

    cdata = smooth3((u-min(min(min(u))))./(max(max(max(u)))-min(min(min(u)))),'box',5);
    fi = smooth3(fi,'box',5);
    p4=patch(isosurface(fi,0));
    isonormals(fi,p4);
    isocolors(cdata,p4);
    set(p4,'FaceColor','interp','EdgeColor','none'),
    camlight, lighting phong
    axis equal, view(-14,40), axis off
    axis([1 Nx 1 Ny 1 Nz]),
    colorbar