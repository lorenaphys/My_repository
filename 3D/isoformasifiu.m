% Program to calculate phase fields in 3 dimensions

%clear all


NF=80;
ep1 = 2;
ep = ep1^2;
sigma = -0.1;
beta = 0.5;
Nx = 40;
Ny = 40;
Nz = 70;
step=200;
dt = 1e-5;
Ab = 0.5;
As = 2;
Af = 2;
Dfi = 1;
Du = 1;
lambda = 1.5;
u1 = 0;
u2 = 1;
u3 = 0;
alpha = 120;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% form of fi and initial values
fi = ones(Nx,Ny,Nz);
r = zeros(Nx,Ny,Nz);

for i=1:Nx
    for j=1:Ny
        for k=1:Nz
      r(i,j,k)=sqrt((i-Nx/2)^2+(j-Ny/2)^2+(k-Nz/2)^2);
      if r(i,j,k)>=15
      fi(i,j,k)=-1;
      end
        end
   
   end
end
%prisma3D
%fiini=fi;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fm = zeros(Nx,Ny,Nz,NF+1);
Um = zeros(Nx,Ny,Nz,NF+1);
Sm = zeros(1,NF+1);

Fm(:,:,:,1) = fi;
fi0=fi;


%%
%Rompimiento de simetria
    R = 9;
    %[~,R1] = min(abs(fi(Nx/2,Ny/2,:)));
    N = 5;
    [X,Y,Z]=meshgrid(1:Nx,1:Ny,1:Nz);
        teta=atan2((Y-Ny/2),(X-Nx/2));
        rad=sqrt((X-Nx/2+.5).^2+(Y-Ny/2+.5).^2);
        %u=1.5*exp(-((X-Nx/2).^2+(Y-Ny/3).^2+(Z-9).^2)/50);
        %u=2.5*rad.*(cos(teta*N)+sin(teta*N)).*(Z/Nz)/max(max(max(rad)))+(exp(-((-X+Nx/2-.5).^2+(-Y+Ny/2-.5).^2+(-Z+R+14).^2)/80));
        u=-2.5*rad.*(cos(teta*N)+sin(teta*N)).*(Z/Nz)/max(max(max(rad)))+(exp(-((X-Nx/2).^2+(Y-Ny/2).^2+(Z-R+14).^2)/50));%pueba 2, N=2
        %u=2.5*rand(Nx,Ny,Nz);
        %u = 1-2*u;
       
     Um(:,:,:,1) = u;
     %Sm(1) = sigma;
    

% [~,R1]=min(abs(fi(Nx/2,Ny/2,:)));

% u = zeros(Nx,Ny,Nz);
% 
%     for i=1:Nx
%         for j=1:Ny
%             for k=1:Nz
%                 u(i,j,k)=1.5*exp(-((i-Nx/3-5)^2+(j-Ny/3-5)^2+(k-13)^2)/20);
%             end
%         end
%     end

u(fi<=-.99)=0;

t = tic(); %comando para visualizar el tiempo de ejecucion
disp(1);

for iter=1:NF
    for iiter=1:step
            
	%[~,bb] = min(abs(fi(Nx/2,Ny/2,:)));
        %R1 = bb-5;
        %u1 = R1;
        %u2 = -R1;

        lapfi=lap3D(fi);
       
        %mu=((fi+ep1*beta.*u).*((fi).^2-1)-ep*lapfi);
        mu = (fi.^2-1).*(fi-ep1*beta*u)-ep*lapfi;
      
        lapmu=lap3D(mu);

	potfi = 2*As*fi.*(fi.^2-1).*((u-u1).^2).*((u-u2).^2) + 2*Af*fi.*((u-u3).^2);
	potu = 2*As*((fi.^2-1).^2).*(u-u1).*(u-u2).*(2*u-u1-u2) + 2*Af*fi.^2.*(u-u3);

        %potffi=2*fi.*(((fi).^2-1)).*((beta*u-R1).^2)+2*fi.*((beta.*u).^2);
        %potfu=2*(((fi).^2-1).^2).*(((beta*u-R1)))+2*(fi.^2).*(beta*u);
        
        lapu=lap3D(u);
        
        F=2*Ab*mu.*(3*fi.^2-1-2*ep1*beta*fi.*u)-2*Ab*ep*lapmu+potfi;
          
        G=-2*Ab*beta*ep1*mu.*(((fi).^2-1))+potu-2*As*lambda*lapu;
         
        Fs=2*((sigma.*lapfi));

        lapF=lap3D(F);
        
        lapFs=lap3D(Fs);

        lapG=lap3D(G);
        
%         B=sum(sum(sum(gradient(fi).*gradient(lapF))));
%         Bs=sum(sum(sum(gradient(fi).*gradient(lapfi))));
%         
%         sigma=B/Bs;
        
        %I=120*(G)*sum(sum(sum((fi>=-.99))))/Nx/Ny/Nz;
        I = alpha*u;
        I(abs(fi)>=.9)=0;
        %I(fi<=0)=0;        
        
        fi=fi+dt*(Dfi*(lapF-lapFs)+I);
        u=u+dt*Du*(lapG);
        %u(fi<=-0.99)=0;
         fi(1,:,:)=fi(2,:,:);
         fi(Nx,:,:)=fi(Nx-1,:,:);
         fi(:,1,:)=fi(:,2,:);
         fi(:,Ny,:)=fi(:,Ny-1,:);         
         fi(:,:,1)=fi(:,:,2);
         fi(:,:,Nz)=fi(:,:,Nz-1);
         %fi(:,:,1)=fi0(:,:,1);
%         %fi(:,:,Nz)=fi(:,:,Nz-1);
%         u(1,:,:)=u(2,:,:);
         u(Nx,:,:)=u(Nx-1,:,:);
         u(:,1,:)=u(:,2,:);
         u(:,Ny,:)=u(:,Ny-1,:);
         u(:,:,1)=u(:,:,2);
         u(:,:,Nz)=u(:,:,Nz-1);
         %noFlux2(fi,u);
        
        %fi(:,:,1)=fiini(:,:,1);

    end
    
    h=isnan(fi(Nx/2,Ny/2,Nz/2));
    if h==1;
        break
    end
	
	Fm(:,:,:,iter+1)=fi;
    Um(:,:,:,iter+1)=u;
    u(fi<=-.99)=0;
    %Sm(iter+1) = sigma;
   % u(fi<=-0.99) = 0;
    disp(iter+1);   
end

time = toc(t)/60;
 

save('oct6h.mat');
