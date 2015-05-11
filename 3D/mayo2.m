% Program to calculate phase fiels in 3 dimensions

clear all

et=1;
dx=1;
NF=100;
sig=0*(1:NF);
ep1=1;
ep=ep1^2;
sigma=1;
Nx=70;
Ny=70;
Nz=50;
R=11;
Afi=1;
N=5;   %this is the cylindrical symmetry



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% form of fi and initial values
%semielipreg;
%formas4cil
%exsemiesf
semiesf3D
%semielip2
%formaneg
%formapinching
%parabola
%cilinder4
%cilindersq2
%formas2
%formas3D
fiini=fi;
% bet=0*fi;
% ra=.5:Nx-.5;
% for i=1:Nx
%     for j=1:Ny
%         rr(i,j)=abs(ra(1,i));
%     end
% end

step=200;
iter=1;
dt=1e-4;
cont=iter;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bet=0*fi;
I=bet;
%%
% [a bb]=min(abs(fi(Nx/2,Ny/2,:)));
% R1=bb;
% 
% for i=1:Nx
%     for j=1:Ny
%         for k=1:Nz
%             bet(i,j,k)=1+3.5*exp(-((i-Nx/2)^2+(j-Ny/2)^2+(k-R1+8)^2)/100);
%             %  bet(i,j,k)=50*k/Nz;
%             %bet(:,:,i)=-3*sin(10*pi*i/Nx-pi)^2+1;
%             % % ee=10;
%         end
%     end
% end
% %bet=14-bet;
% bet(fi<=-.9);
% %load betu
% %bet=10*u;
% %bet(1:Nx,1:Ny)=3000*u(1:Nx,1:Ny);
% 
% bet(fi<=-.99)=0;
% figure(3)
% clf
% cdata =
% smooth3((bet-min(min(min(bet))))./(max(max(max(bet)))-min(min(min(bet)))),'box',5);
% %
% %      p3=patch(isosurface(bet,0));
% %      isonormals(bet,p3);
% %      isocolors(cdata,p3);
% %      hold
% [x,y,z] = meshgrid(1:1:Nx,1:1:Ny,1:1:Nz);
% xslice = [Nx/2,Ny/2];yslice = Ny/2; zslice = [0,10];
% p3=slice(x,y,z,bet,xslice,yslice,zslice);
% set(p3,'FaceColor','interp','EdgeColor','none'),
% axis equal, view(-60,8)
% pause(.01)
     fi = smooth3(fi,'box',5);
fix0(:,:)=fi(Nx/2,:,:);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

%%
%load phase3d-2
for iter=cont:NF
    for iiter=1:step
        
[a bb]=min(abs(fi(Nx/2,Ny/2,:)));
        R1=bb;
        
%         for i=1:Nx
%             for j=1:Ny
%                 for k=1:Nz
% %                     bet1(i,j,k)=exp(-((i-Nx/2+9)^2+(j-Ny/2)^2+(k-(R1-4))^2)/10);
% %                     bet2(i,j,k)=exp(-((i-Nx/2-9)^2+(j-Ny/2)^2+(k-(R1-4))^2)/10);
%                     %bet(i,j,k)=1+3.5*exp(-((i-Nx/2)^2+(j-Ny/2)^2+(k-(R1-8)^2))/90);
%                       %bet(i,j,k)=0;%R*(exp(k/R1-1)-exp(1/R1-1));
%                     %bet(:,:,i)=-3*sin(10*pi*i/Nx-pi)^2+1;
%                     % % ee=10;
%                 end
%             end
%         end
        
        [X,Y,Z]=meshgrid(1:Nx,1:Ny,1:Nz);
        teta=atan2((Y-Ny/2),(X-Nx/2));
        rad=sqrt((X-Nx/2+.5).^2+(Y-Ny/2+.5).^2);
        bet=2.5*rad.*(cos(teta*N)+sin(teta*N)).*(Z/Nz)/max(max(max(rad)))+exp(-((X-Nx/2).^2+(Y-Ny/2).^2+(Z-(R1)).^2)/60);
        bet=3*(bet);
        %bet=3.5*(bet-min(min(min(bet))))+1;
        %bet=14-bet;
        %bet(fi<=-.9);
        %load betu
        %bet=10*u;
        %bet(1:Nx,1:Ny)=3000*u(1:Nx,1:Ny);
        %bet=3.5*(bet1+bet2)+.6;
        
                bet(fi<=-ep)=0;
                bet(fi>=ep)=0;
%%                
        H=fi;
        lap3D
        lapfi=lapH;
        
        
        
        mu=((fi-ep1.*bet).*((fi).^2-1)-ep*lapfi);
        
        H=mu;
        lap3D
        lapmu=lapH;
        
        
        F=Afi*((3*fi.^2-1-2*ep1*fi.*bet).*mu-ep*lapmu);
        
        
        
        H=F;
        lap3D
        lapF=lapH;
        
        
        Fs=-sigma*lapfi;
        
        
        H=Fs;
        lap3D
        lapFs=lapH;
        
        
       %I=120*(lapF +lapFs);
       
        %       fi=fi-dt*(F+Fs);
        I=100*(F)*sum(sum(sum((fi>=-.99))))/Nx/Ny/Nz;
        I(fi<=0)=0;        
        fi=fi+dt*(lapF+lapFs-I);
%        fi(:,:,1)=fiini(:,:,1);

    end
    
    %%
    sig(iter)=sigma;
    
    h=isnan(fi(Nx/2,Ny/2,Nz/2));
    if h==1;
        break
    end
    sigma;
    
    
    
  %iter
  
  fi = smooth3(fi,'box',5);
  bet = smooth3(bet,'box',5);
  
  
  Fm(:,:,:,iter)=fi(:,:,:);
  U(:,:,iter)=bet(:,:);
  Um(:,:,:,iter)=bet(:,:,:); 
    %%
            %%
iter            
figure(1)
clf
isosurface(fi,0)
view(64,26),
axis equal
axis([1 Nx 1 Ny 1 Nz])
%%
%% 
      figure(2)
%     clf
%     cdata =
%     smooth3((bet-min(min(min(bet))))./(max(max(max(bet)))-min(min(min(bet
%     )))),'box',5);
   

%cdata = smooth3((bet-min(min(min(bet))))./(max(max(max(bet)))-min(min(min(bet)))),'box',5);
[x,y,z] = meshgrid(1:1:Nx,1:1:Ny,1:1:Nz);
xslice = [Nx/2,Nx/2];yslice = [Ny/2-10,Ny/2+10]; zslice = [R,2*R];
p3=slice(x,y,z,bet,xslice,yslice,zslice);
set(p3,'FaceColor','interp','EdgeColor','none','FaceAlpha',0.5);
axis equal, view(68,16), 
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
   axis equal
%    getframe(gcf);
   hold off
    
end
save may7-1
%exit
