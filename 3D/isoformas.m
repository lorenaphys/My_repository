% Program to calculate phase fiels in 3 dimensions

clear all

et=1;
NF=200;
ep1=1;
ep=ep1^2;
sigma=1;
k=0.5;
kb=3*k*(sqrt(2))/4/ep/ep1;
Nx=80;
Ny=50;
Nz=50;
R=11;
Afi=1;
%N=5;   %this is the cylindrical symmetry
v = 0.05; %reduce volume
V = 4000; %Prisma volume



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% form of fi and initial values

prisma3D

fiini=fi;


step=200;
dt=1e-4;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%bet=0*fi;
%I=bet;

%fix0(:,:)=fi(:,:,Nz/2);

%%
ao = 4*pi*(3*V/(4*pi*v))^(2/3);
%vol=.46; 

%ao=((1-fi.*fi).^2);

Fm = zeros(Nx,Ny,Nz,NF);


for iter=1:NF
    for iiter=1:step

		lapfi = lap3D(fi);

        mu=((fi).*((fi).^2-1)-ep*lapfi);
        
        lapmu = lap3D(mu);

        F=Afi*((3*fi.^2-1).*mu-ep*lapmu);
 
        gradfi=gradient(fi);
        
        [fix, fiy, fiz] = gradient(fi);

        %sigma= sigma+dt*((3/(4*sqrt(2)*ep1))*((1-fi.*fi).^2)-ao);
        
        sigma= sigma+dt*((3/(4*sqrt(2)*ep1))*(fix.^2 + fiy.^2 + fiz.^2)-ao);
        
        Fs=((3*sqrt(2)/(4*ep1*kb))*(sigma.*lapfi+gradient(sigma).*gradfi));

		lapF = lap3D(F);
		
		lapFs = lap3D(Fs);

        %I=gradient(F+Fs).^2;

        fi=fi+dt*(lapF+lapFs);

        fi(:,:,1)=fiini(:,:,1);
        
        %fi(1,:,:)=fi(2,:,:);
        %fi(Nx,:,:)=fi(Nx-1,:,:);
        %fi(:,1,:)=fi(:,2,:);
        %fi(:,Ny,:)=fi(:,Ny-1,:);
        %fi(:,:,1)=fi(:,:,2);
        %fi(:,:,Nz)=fi(:,:,Nz-1);

    end

%     h=isnan(fi(Nx/2,Ny/2,Nz/2));
%     if h==1;
%         break
%     end
	Fm(:,:,:,iter)=fi(:,:,:);
	%save('test.txt','iter')
 
end

save agosto11a;

exit

