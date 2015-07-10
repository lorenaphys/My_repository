% Program to calculate phase fields in 3 dimensions

clear all

NF=200;
%sig=0*(1:NF);
ep1=2;
ep=ep1^2;
sigma=0.1;
beta=0.5;
Nx=40;
Ny=40;
Nz=70;
step=200;
dt=1e-5;
Ab = 0.5;
As = 2;
Af = 2;
Dfi = 1;
Du = 1;
lambda = 0.1;
u1 = 0;
u2 = 1;
u3 = 0;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% form of fi and initial values
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
%prisma3D
%fiini=fi;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%u = zeros(Nx,Ny,Nz);
%I = zeros(Nx,Ny,Nz);
Fm = zeros(Nx,Ny,Nz,NF);
U = zeros(Nx,Ny,Nz,NF);

%%
%Rompimiento de simetria
R = 9;
    %N = 6;
    [X,Y,Z]=meshgrid(1:Nx,1:Ny,1:Nz);
        %teta=atan2((Y-Ny/2),(X-Nx/2));
        %rad=sqrt((X-Nx/2+.5).^2+(Y-Ny/2+.5).^2);
        u=1.5*exp(-((X-Nx/3-.5).^2+(Y-Ny/3-5-.5).^2+(Z-R+2).^2)/50);
        %u=2.5*rad.*(cos(teta*N)+sin(teta*N)).*(Z/Nz)/max(max(max(rad)))+(exp(-((-X+Nx/2-.5).^2+(-Y+Ny/2-.5).^2+(-Z+R+14).^2)/80));
        %u=-2.5*rad.*(cos(teta*N)+sin(teta*N)).*(Z/Nz)/max(max(max(rad)))+(exp(-((X-Nx/2).^2+(Y-Ny/2).^2+(Z-R+14).^2)/50));%pueba 2, N=2
        %u=2.5*rand(Nx,Ny,Nz);
    

% [~,R1]=min(abs(fi(Nx/2,Ny/2,:)));
% 
% for i=1:Nx
%     for j=1:Ny
%         for k=1:Nz
%             u(i,j,k)=1+3.5*exp(-((i-Nx/2)^2+(j-Ny/2)^2+(k-R1+8)^2)/200);%primera prueba ancho = 200;
%         end
%     end
% end

%u(fi<=-.9);

u(fi<=-.99)=0;

for iter=1:NF
    for iiter=1:step
            
		%[~,bb]=min(abs(fi(Nx/2,Ny/2,:)));
         %R1=bb-5;

        lapfi=lap3D(fi);
       
        mu=((fi+ep1*beta.*u).*((fi).^2-1)-ep*lapfi);
      
        lapmu=lap3D(mu);

		potffi = 2*As*fi.*(fi.^2-1).*((u-u1).^2).*((u-u2).^2) + 2*Af*fi.*((u-u3).^2);
		potfu = 2*As*((fi.^2-1).^2).*(u-u1).*(u-u2).*(2*u-u1-u2) + 2*Af*fi.^2.*(u-u3);

        %potffi=2*fi.*(((fi).^2-1)).*((beta*u-R1).^2)+2*fi.*((beta.*u).^2);
        %potfu=2*(((fi).^2-1).^2).*(((beta*u-R1)))+2*(fi.^2).*(beta*u);
        
        lapu=lap3D(u);
        
        F=2*Ab*mu.*(3*fi.^2-1-2*ep1*beta*fi.*u)-2*Ab*ep*lapmu+potffi;
          
        G=-2*Ab*beta*ep1*mu.*(((fi).^2-1))+potfu-2*As*lambda*lapu;
         
        Fs=2*((sigma.*lapfi));

        lapF=lap3D(F);
        
        lapFs=lap3D(Fs);

        lapG=lap3D(G);
        
        B=sum(sum(sum(gradient(fi).*gradient(lapF))));
        Bs=sum(sum(sum(gradient(fi).*gradient(lapfi))));
        
        sigma=B/Bs;
        %sig(iter)=sigma;
        
        I=120*(G)*sum(sum(sum((fi>=-.99))))/Nx/Ny/Nz;
        I(fi<=0)=0;        
        
        fi=fi+dt*(Dfi*(lapF-lapFs)+I);
        u=u+dt*Du*(lapG);
        
        fi(1,:,:)=fi(2,:,:);
        fi(Nx,:,:)=fi(Nx-1,:,:);
        fi(:,1,:)=fi(:,2,:);
        fi(:,Ny,:)=fi(:,Ny-1,:);
        %fi(:,:,1)=fi(:,:,2);
        %fi(:,:,Nz)=fi(:,:,Nz-1);
        u(1,:,:)=u(2,:,:);
        u(Nx,:,:)=u(Nx-1,:,:);
        u(:,1,:)=u(:,2,:);
        u(:,Ny,:)=u(:,Ny-1,:);
        %u(:,:,1)=u(:,:,2);
        %u(:,:,Nz)=u(:,:,Nz-1);
        
        
        %fi(:,:,1)=fiini(:,:,1);

    end
    
    h=isnan(fi(Nx/2,Ny/2,Nz/2));
    if h==1;
        break
    end

	Fm(:,:,:,iter)=fi(:,:,:);
    U(:,:,:,iter)=u(:,:,:);
    save('test','iter')
	%disp(iter)
   
end

	save julio9j;

exit
