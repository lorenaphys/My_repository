% Program to calculate phase fiels in 3 dimensions

clear all

et=1;
dx=1;
NF=150;
%sig=0*(1:NF);
ep1=2;
ep=ep1^2;
sigma=-0.1;
beta=0.5;
gamma=0.0;
D=1.5;
Dfi=1;
Du=30;
Nx=60;
Ny=60;
Nz=60;
R=11;
Afi=0.25;
As=1;
Af=1;
N=5;   %this is the cylindrical symmetry
step=300;
iter=1;
dt=1e-5;
cont=iter;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% form of fi and initial values
semiesf3D
%prisma3D
fiini=fi;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u=0*fi;
I=u;
potffi=u;
potfu=u;
str=u;
Ft=u;

%   No Gaussicana
[a bb]=min(abs(fi(Nx/2,Ny/2,:)));
R1=bb;

for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            u(i,j,k)=1+3.5*exp(-((i-Nx/2)^2+(j-Ny/2)^2+(k-R1+8)^2)/20);
            %  bet(i,j,k)=50*k/Nz;
            %bet(:,:,i)=-3*sin(10*pi*i/Nx-pi)^2+1;
            % % ee=10;
        end
    end
end

%u=0.1*ones(Nx,Ny,Nz);
%u=ones(Nx,Ny,Nz);
%bet=14-bet;
%u(fi<=-.9);
%load betu
%bet=10*u;
%bet(1:Nx,1:Ny)=3000*u(1:Nx,1:Ny);

%u(fi<=-.99)=0;
%figure(3)
%clf
%cdata = smooth3((u-min(min(min(u))))./(max(max(max(u)))-min(min(min(u)))),'box',5);
%
%      p3=patch(isosurface(bet,0));
%      isonormals(bet,p3);
%      isocolors(cdata,p3);
%      hold
%[x,y,z] = meshgrid(1:1:Ny,1:1:Nx,1:1:Nz);
%xslice = [Nx/2,Nx/2];yslice = [Ny/2,Ny/2]; zslice = [0,1];
%xslice = [Nx/2,Ny/2];yslice = Ny/2; zslice = [0,10];
%p3=slice(x,y,z,u,xslice,yslice,zslice);
%set(p3,'FaceColor','interp','EdgeColor','none'),
%axis equal, view(56,30)
%pause(.01)
%fi = smooth3(fi,'box',5);
fix0(:,:)=fi(:,Ny/2,:);
      
%%  Gaussiana
% [a bb]=min(abs(fi(Nx/2,Ny/2,:)));
% R1=bb;
% 
% for i=1:Nx
%     for j=1:Ny
%         for k=1:Nz
%             u(i,j,k)=1+3.5*exp(-((i-Nx/2)^2+(j-Ny/2)^2+(k-R1+8)^2)/100);
%             %  bet(i,j,k)=50*k/Nz;
%             %bet(:,:,i)=-3*sin(10*pi*i/Nx-pi)^2+1;
%             % % ee=10;
%         end
%     end
% end
% %bet=14-bet;
% u(fi<=-.9);
% %load betu
% %bet=10*u;
% %bet(1:Nx,1:Ny)=3000*u(1:Nx,1:Ny);
% 
% u(fi<=-.99)=0;
% figure(3)
% clf
% cdata = smooth3((u-min(min(min(u))))./(max(max(max(u)))-min(min(min(u)))),'box',5);
% %
% %      p3=patch(isosurface(bet,0));
% %      isonormals(bet,p3);
% %      isocolors(cdata,p3);
% %      hold
% [x,y,z] = meshgrid(1:1:Nx,1:1:Ny,1:1:Nz);
% xslice = [Nx/2,Ny/2];yslice = Ny/2; zslice = [0,10];
% p3=slice(x,y,z,u,xslice,yslice,zslice);
% set(p3,'FaceColor','interp','EdgeColor','none'),
% axis equal, view(-60,8)
% pause(.01)
% fi = smooth3(fi,'box',5);
% fix0(:,:)=fi(Nx/2,:,:);





%%
%load phase3d-2
for iter=cont:NF
    for iiter=1:step
%         
%  [a bb]=min(abs(fi(Nx/2,Ny/2,:)));
%          %R1=bb-5;
%          R1=bb;
 %%        
        
%         [X,Y,Z]=meshgrid(1:Nx,1:Ny,1:Nz);
%         teta=atan2((Y-Ny/2),(X-Nx/2));
%         rad=sqrt((X-Nx/2+.5).^2+(Y-Ny/2+.5).^2);
%         u=2.5*rad.*(cos(teta*N)+sin(teta*N)).*(Z/Nz)/max(max(max(rad)))+exp(-((X-Nx/2).^2+(Y-Ny/2).^2+(Z-(R1)).^2)/40);
%         u=3*(u);
% 
%         
%                 u(fi<=-.5)=0;
%                 u(fi>=.5)=0;
                %%

        lapfi=lapf3D(fi);
       
        mu=Afi*((fi+ep1*beta*u).*((fi).^2-1)-ep*lapfi);
      
        lapmu=lapf3D(mu);

         potffi=As*(2*fi.*(((fi).^2-1)).*((beta*u-0).^2).*((beta.*u-1).^2))+Af*2*fi.*((beta.*u-0).^2);
         potfu=As*(2*(((fi).^2-1).^2).*(((beta*u-0).^2).*(beta*u-1)+(beta*u-0).*((beta*u-1).^2)))+Af*2*(fi.^2).*(beta*u-0);
%         potffi=2*fi.*(((fi).^2-1)).*((beta*u-R1).^2)+2*fi.*((beta.*u-0).^2);
%         potfu=2*(((fi).^2-1).^2).*(((beta*u-R1)))+2*(fi.^2).*(beta*u-0);
        
        
        
        lapu=lapf3D(u);
        
        F=Afi*(2*(3*fi.^2-1-2*ep1*beta*fi.*u).*mu-2*ep*lapmu)+potffi-gamma*lapu;
          
        G=Afi*(2*ep1*mu.*(((fi).^2-1)))+potfu-gamma*lapfi-D*lapu;
         
        Fs=1*((sigma.*lapfi));

        lapF=lapf3D(F);
        
        lapFs=lapf3D(Fs);

        lapG=lapf3D(G);
        
        %I=sum(sum(sum(gradient(fi).*gradient(lapF))));
        %Is=sum(sum(sum(gradient(fi).*gradient(lapfi))));
        
        %sigma=I/Is;
        %sig(iter)=sigma;
        
        %If=20*(F)*sum(sum(sum((fi>=-.99))))/Nx/Ny/Nz;
        I=20*sum(sum(sum(F.*fi)))/Nx/Ny/Nz;
        %If(fi<=0)=0; 
        
        Ft=-fi.*gradient(mu);
         
        fi=fi+dt*Dfi*(lapF+lapFs)+dt*I*(fi>=-.9);
        u=u+dt*Du*(lapG);
        %u=u+dt*gradient(str.*gradient(G));
%         fi(1,:,:)=fi(2,:,:);
%         fi(Nx,:,:)=fi(Nx-1,:,:);
%         fi(:,1,:)=fi(:,2,:);
%         fi(:,Ny,:)=fi(:,Ny-1,:);
         fi(:,:,1)=fi(:,:,2);
%         fi(:,:,Nz)=fi(:,:,Nz-1);
%         u(1,:,:)=u(2,:,:);
%         u(Nx,:,:)=u(Nx-1,:,:);
%         u(:,1,:)=u(:,2,:);
%         u(:,Ny,:)=u(:,Ny-1,:);
%         u(:,:,1)=u(:,:,2);
%         u(:,:,Nz)=u(:,:,Nz-1);


%	fi(:,:,1)=0; 
        
%        fi(:,:,1)=fiini(:,:,1);

        %u(fi<=-.99)=0;
    end
    
    %%
    str1o=mu.^2+0.5*sigma*gradient(fi).*gradient(fi)+0.5*gradient(u).*gradient(u)-fi.*F;
    str=zeros(Nx,Ny,Nz);
    for hdo=1:Nx
    str1=str1o(hdo,hdo,hdo);
    end
    
    str2=-2*sigma*gradient(fi).*gradient(fi);
    str3=-2*ep*gradient(mu).*gradient(fi);
    str4=-2*ep*mu.*gradient(fi).*gradient(fi);
    
    str=str1+str2+str3+str4;
    
    %%
    %sig(iter)=sigma;
    
    h=isnan(fi(Nx/2,Ny/2,Nz/2));
    if h==1;
        break
    end
    sigma;
    
    %sum(sum(sum((fi>=0))))
    
  iter
  
  fi = smooth3(fi,'box',5);
  %sigma = smooth3(sigma,'box',5);
  u = smooth3(u,'box',5);
  
  
  Fm(:,:,:,iter)=fi(:,:,:);
  %U(:,:,iter)=bet(:,:);
  Um(:,:,:,iter)=u(:,:,:); 
  Strm(:,:,:,iter)=str(:,:,:);
  Ftm(:,:,:,iter)=Ft(:,:,:); 
  %%
            %%
%iter            
figure(1)
clf
isosurface(fi,0)
view(45,10),
axis equal
axis([1 Ny 1 Nx 1 Nz])
%%
%% 
%      figure(2)
%     clf
%%     cdata =
%%     smooth3((bet-min(min(min(bet))))./(max(max(max(bet)))-min(min(min(bet
%%     )))),'box',5);
   

%%cdata = smooth3((bet-min(min(min(bet))))./(max(max(max(bet)))-min(min(min(bet)))),'box',5);
[x,y,z] = meshgrid(1:1:Ny,1:1:Nx,1:1:Nz);
xslice = [Nx/2,Nx/2];yslice = [Ny/2,Ny/2]; zslice = [0,1];
p3=slice(x,y,z,fi,xslice,yslice,zslice);
set(p3,'FaceColor','interp','EdgeColor','none','FaceAlpha',0.5);
axis equal, view(56,30), 
colorbar
%%pause(.01)
% 
%           M(:,:,iter)=getframe;

%para hacer peliculas 3D usar cine3D


    %%
  figure(3)
  clf
  fix(:,:)=fi(:,Ny/2,:);
  contour(fix,[0 0],'k')
  hold on
  contour(fix0,[0 0],'r')
  view(-90,90)
  axis equal
%%    getframe(gcf);
%   hold off
   
   %%
        figure(4)
    clf
xslice = [25,20];yslice = [Ny/2,Ny/2]; zslice = [0,1];
p3=slice(x,y,z,str,xslice,yslice,zslice);
set(p3,'FaceColor','interp','EdgeColor','none','FaceAlpha',0.5);
axis equal, view(-100,30), 
colorbar

%%
figure(5)
clf
   fiss=fi;  
   cdata = smooth3((u-min(min(min(u))))./(max(max(max(u)))-min(min(min(u)))),'box',5);
   fiss = smooth3(fiss,'box',5);
   p4=patch(isosurface(fiss,0));
   isonormals(fiss,p4);
   isocolors(cdata,p4);
   set(p4,'FaceColor','interp','EdgeColor','none'),
   camlight, lighting phong
   axis equal, view(-14,40), axis off
   axis([1 Nx 1 Ny 1 Nz/2]),
   colorbar
    
end
save isoformas-may16
%exit
