% Program to calculate phase fields in 3 dimensions

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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% form of fi and initial values

%prisma3D
esf3D

a = area(fi,ep1);
ao = a;

fiini=fi;


step=100;
dt=1e-4;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v = reducido(fi,ep1);

Fm = zeros(Nx,Ny,Nz,NF+1);
Fm(:,:,:,1) = fi;
Vm = zeros(1,NF+1);
Vm(1) = v;


for iter=1:NF
    for iiter=1:step

		lapfi = lap3D(fi);

        mu=((fi).*((fi).^2-1)-ep*lapfi);
        
        lapmu = lap3D(mu);

        F=Afi*((3*fi.^2-1).*mu-ep*lapmu);
 
        gradfi=gradient(fi);
        
        a = area(fi,ep1);
        
        sigma= sigma+dt*(a-ao);
        
        Fs=((3*sqrt(2)/(4*ep1*kb))*(sigma.*lapfi+gradient(sigma).*gradfi));

		lapF = lap3D(F);
		
		lapFs = lap3D(Fs);

        %I=gradient(F+Fs).^2;
        
        v = reducido(fi,ep1);

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
	Fm(:,:,:,iter+1)=fi(:,:,:);
    Vm(iter+1) = v;
	%save('test.txt','iter')
 
end

save('oct21b.mat','Fm','Vm');


