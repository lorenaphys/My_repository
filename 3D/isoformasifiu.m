% Program to calculate phase fields in 3 dimensions

clear all

NF=5000;
%sig=0*(1:NF);
ep1=.1;
ep=ep1^2;
sigma=-5.5;
beta=0.5;
gamma=0;
Nx=60;
Ny=40;
Nz=40;
R=11;
Afi=1;
N=5;   %this is the cylindrical symmetry
step=200;
iter=1;
dt=1e-4;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% form of fi and initial values
fi=ones(Nx,Ny,Nz);

for i=1:Nx
    for j=1:Ny
        for k=1:Nz
      r(i,j,k)=sqrt((i-Nx/2)^2+(j-Ny/2)^2+(k)^2);
      if r(i,j,k)>=R
      fi(i,j,k)=-1;
      end
        end
   
   end
end
%prisma3D
fiini=fi;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = zeros(Nx,Ny,Nz);
I = u;

%%
[~,R1]=min(abs(fi(Nx/2,Ny/2,:)));

for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            u(i,j,k)=1+3.5*exp(-((i-Nx/2)^2+(j-Ny/2)^2+(k-R1+8)^2)/100);
        end
    end
end

%u(fi<=-.9);

u(fi<=-.99)=0;

for iter=1:NF
    for iiter=1:step
            
		[~,bb]=min(abs(fi(Nx/2,Ny/2,:)));
         R1=bb-5;

        lapfi=lap3D(fi);
       
        mu=((fi+ep1*beta.*u).*((fi).^2-1)-ep*lapfi);
      
        lapmu=lap3D(mu);

        potffi=2*fi.*(((fi).^2-1)).*((beta*u-R1).^2)+2*fi.*((beta.*u).^2);
        potfu=2*(((fi).^2-1).^2).*(((beta*u-R1)))+2*(fi.^2).*(beta*u);
        
        lapu=lap3D(u);
        
        F=2*(3*fi.^2-1-2*ep1*beta*fi.*u).*mu-2*ep*lapmu+potffi-gamma*beta*lapu;
          
        G=2*mu*ep1.*(((fi).^2-1))+potfu-gamma*lapfi;
         
        Fs=1*((sigma.*lapfi));

        lapF=lap3D(F);
        
        lapFs=lap3D(Fs);

        lapG=lap3D(G);
        
        I=sum(sum(sum(gradient(fi).*gradient(lapF))));
        Is=sum(sum(sum(gradient(fi).*gradient(lapfi))));
        
        sigma=I/Is;
        %sig(iter)=sigma;
        
        I=120*(F)*sum(sum(sum((fi>=-.99))))/Nx/Ny/Nz;
        I(fi<=0)=0;        
        
        fi=fi+dt*(lapF-lapFs-I);
        u=u+dt*(lapG);
        
        fi(1,:,:)=fi(2,:,:);
        fi(Nx,:,:)=fi(Nx-1,:,:);
        fi(:,1,:)=fi(:,2,:);
        fi(:,Ny,:)=fi(:,Ny-1,:);
        fi(:,:,1)=fi(:,:,2);
        fi(:,:,Nz)=fi(:,:,Nz-1);
        u(1,:,:)=u(2,:,:);
        u(Nx,:,:)=u(Nx-1,:,:);
        u(:,1,:)=u(:,2,:);
        u(:,Ny,:)=u(:,Ny-1,:);
        u(:,:,1)=u(:,:,2);
        u(:,:,Nz)=u(:,:,Nz-1);
        
        
        fi(:,:,1)=fiini(:,:,1);

    end
    
    h=isnan(fi(Nx/2,Ny/2,Nz/2));
    if h==1;
        break
    end
    
	disp(iter)
   
end

	save(['/iter' num2str(iter)]);

exit