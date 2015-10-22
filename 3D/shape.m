% Program to calculate phase fields in 3 dimensions

%clear all

%et=1;
NF=200;
ep1=1;
ep=ep1^2;
sigma=0.01;
k=0.5;
kb=3*k*(sqrt(2))/4/ep/ep1;
Nx=80;
Ny=50;
Nz=50;
%R=11;
Afi=1;
%N=5;   %this is the cylindrical symmetry 
V = 4000; %Prisma volume



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% form of fi and initial values

prisma3D

a = area(fi);
ao = a;

fiini=fi;


step=100;
dt=1e-4;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%bet=0*fi;
%I=bet;

%fix0(:,:)=fi(:,:,Nz/2);

%%
%ao = 4*pi*(3*V/(4*pi*v))^(2/3);
Ro = sqrt(a/4*pi);
v = V/(4*pi/3)*Ro^3;
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
        
        a = area(fi);
        
        sigma= sigma+dt*(a-ao);
        
        V = sum(sum(sum(fi>0.99)));
        
        Ro = sqrt(a/4*pi);
        v = V/(4*pi/3)*Ro^3;
        
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

save('oct16a.mat','Fm','v');


