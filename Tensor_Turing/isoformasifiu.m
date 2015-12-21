% Program to calculate phase fiels in 3 dimensions

clear all

et=1;
dx=1;
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
cont=iter;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% form of fi and initial values
semiesf3D
%prisma3D
fiini=fi;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u=0*fi;
I=u;

%%
[a bb]=min(abs(fi(Nx/2,Ny/2,:)));
R1=bb;

for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            u(i,j,k)=1+3.5*exp(-((i-Nx/2)^2+(j-Ny/2)^2+(k-R1+8)^2)/100);
            %  bet(i,j,k)=50*k/Nz;
            %bet(:,:,i)=-3*sin(10*pi*i/Nx-pi)^2+1;
            % % ee=10;
        end
    end
end
%bet=14-bet;
u(fi<=-.9);
%load betu
%bet=10*u;
%bet(1:Nx,1:Ny)=3000*u(1:Nx,1:Ny);

u(fi<=-.99)=0;
figure(3)
clf
%cdata = smooth3((u-min(min(min(u))))./(max(max(max(u)))-min(min(min(u)))),'box',5);
%
%      p3=patch(isosurface(bet,0));
%      isonormals(bet,p3);
%      isocolors(cdata,p3);
%      hold
[x,y,z] = meshgrid(1:1:Ny,1:1:Nx,1:1:Nz);
xslice = [Nx/2-15,Nx/2-15];yslice = [Ny/2+15,Ny/2+15]; zslice = [0,1];
%xslice = [Nx/2,Ny/2];yslice = Ny/2; zslice = [0,10];
p3=slice(x,y,z,u,xslice,yslice,zslice);
set(p3,'FaceColor','interp','EdgeColor','none'),
axis equal, view(56,30)
%pause(.01)
%fi = smooth3(fi,'box',5);
fix0(:,:)=fi(:,Ny/2,:);
      
%%
%load phase3d-2
for iter=cont:NF
    for iiter=1:step
%         
 [a bb]=min(abs(fi(Nx/2,Ny/2,:)));
         R1=bb-5;

        lapfi=lapf3D(fi);
       
        mu=((fi+ep1*beta.*u).*((fi).^2-1)-ep*lapfi);
      
        lapmu=lapf3D(mu);

%         potffi=2*fi.*(((fi).^2-1)).*((beta*u-R1).^2).*((beta.*u+R1).^2)+2*fi.*((beta.*u-0).^2);
%         potfu=2*(((fi).^2-1).^2).*(((beta*u-R1).^2).*(beta*u+R1)+(beta*u-R1).*((beta*u+R1).^2))+2*(fi.^2).*(beta*u-0);
        potffi=2*fi.*(((fi).^2-1)).*((beta*u-R1).^2)+2*fi.*((beta.*u-0).^2);
        potfu=2*(((fi).^2-1).^2).*(((beta*u-R1)))+2*(fi.^2).*(beta*u-0);
        
        lapu=lapf3D(u);
        
        F=2*(3*fi.^2-1-2*ep1*beta*fi.*u).*mu-2*ep*lapmu+potffi-gamma*beta*lapu;
          
        G=2*mu*ep1.*(((fi).^2-1))+potfu-gamma*lapfi;
         
        Fs=1*((sigma.*lapfi));

        lapF=lapf3D(F);
        
        lapFs=lapf3D(Fs);

        lapG=lapf3D(G);
        
        I=sum(sum(sum(gradient(fi).*gradient(lapF))));
        Is=sum(sum(sum(gradient(fi).*gradient(lapfi))));
        
        sigma=I/Is;
        sig(iter)=sigma;
        
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
    
    %%
    %sig(iter)=sigma;
    
    h=isnan(fi(Nx/2,Ny/2,Nz/2));
    if h==1;
        break
    end
    sigma;
    
    sum(sum(sum((fi>=0))))
    
  %iter
  
  fi = smooth3(fi,'box',5);
  %sigma = smooth3(sigma,'box',5);
  u = smooth3(u,'box',5);
  
  
  Fm(:,:,:,iter)=fi(:,:,:);
  %U(:,:,iter)=bet(:,:);
  Um(:,:,:,iter)=u(:,:,:); 
    %%
            %%
iter            
figure(1)
clf
isosurface(fi,0)
view(45,10),
axis equal
axis([1 Ny 1 Nx 1 Nz])
%%
%% 
      figure(2)
     clf
%     cdata =
%     smooth3((bet-min(min(min(bet))))./(max(max(max(bet)))-min(min(min(bet
%     )))),'box',5);
   

%cdata = smooth3((bet-min(min(min(bet))))./(max(max(max(bet)))-min(min(min(bet)))),'box',5);
[x,y,z] = meshgrid(1:1:Ny,1:1:Nx,1:1:Nz);
xslice = [Nx/2-15,Nx/2-15];yslice = [Ny/2+15,Ny/2+15]; zslice = [0,1];
p3=slice(x,y,z,u,xslice,yslice,zslice);
set(p3,'FaceColor','interp','EdgeColor','none','FaceAlpha',0.5);
axis equal, view(56,30), 
colorbar
%pause(.01)
% 
%           M(:,:,iter)=getframe;

%para hacer peliculas 3D usar cine3D


    %%
   figure(3)
   fix(:,:)=fi(:,Ny/2,:);
   contour(fix,[0 0],'k')
   hold on
   contour(fix0,[0 0],'r')
   view(-90,90)
   axis equal
%    getframe(gcf);
   hold off
    
end
save may8
%exit
