% Program to calculate phase fiels in 3 dimensions

%clear all

dx=1;
NF=300;
sig=0*(1:NF);
ep1=1;
ep=ep1^2;
sigma=0.1;
Nx=40;
Ny=40;
Nz=70;
R=11;
N=3;
alpha = 50;
%%% wfi=0.5*Afi*mu.^2 + 0.5*sigma*gradfi + 0.5*AI*UI + 0.5*As*BI + 0.5*As*(fi.^2 -1).*BS + 0.5*Af*(fi.^2).*BF;
%%Energy=chemical+ surface tension+ fi-u-interaction + u-s interaction +
%%membrane interaction near + far
%%%%%%%%%%%%%%%%%%Difusion constans %%%%%%%%%%%%%%%%%%%%%%555
Du1=10;
Du2=0.0;
Du3=0.00;
Du4=0;
Dfi=1;
eta=5;
%%%%%%%%%%%%%%%%%%%%% strength of the fields  %%%%%%%%%%%
As=0; 
AI=0; 
Afi=1;
Af=0;
%%%%%%%%%%%%%%%%%% spontaneous interaction constans %%%%%%%%%%%%%%
bet1=.5;
bet2=0;
bet3=0;
bet4=0;
%%%%%%%%%%%%%%%%%%%% fi and u's interaction constans  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gam1=0;
gam2=0;
gam3=0;
gam4=0;
%%%%%%%%%%%%%%%%%%%%%%% u's surface tension constants   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d11=0.07;
d12=0;
d13=0;
d14=0;
d21=0;
d22=0;
d23=0;
d24=0;
d31=0;
d32=0;
d33=0;
d34=0;
d41=0;
d42=0;
d43=0;
d44=0;
%%%%%%%%%%%%%%%  fixed points %%%%%%%%%%%%%%%
u11=0;%-0.45;
u12=0;%-2.72;
u13=0;%-4.07;
u14=0;%-2.36;
u21=1;%-1.35;
u22=0;%-2.32;
u23=0;%1.01;
u24=0;%1.82;
u31=0;%-0.5;
u32=0;%-0.5;
u33=0;%0.5;
u34=0;%0.5;
u41=0;%-0.5;
u42=0;%0.5;
u43=0;%-0.5;
u44=0;%0.5;

u1f=0;
u2f=0;
u3f=0;
u4f=0;
%load dominio
semiesf3D
%formas3D
%prisma3D
fiini=fi;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initial conditions for u and beta %%%%%

% u1=1; 
% u2=-.0; 
% 
% uout=0.0; 



        
[X,Y,Z]=meshgrid(1:Nx,1:Ny,1:Nz);
%         
        tetar=0;   % X rotation
        fir=0;     % Z rotation
        RX=(X)*cos(fir)-(Y)*sin(fir)*cos(tetar)+(Z)*sin(fir)*sin(tetar);
        RY=(X)*sin(fir)+(Y)*cos(fir)*cos(tetar)-(Z)*cos(fir)*sin(tetar);
        RZ=(Y)*sin(tetar)+(Z)*cos(tetar);
        teta=atan2((RY-Ny/2),(RX-Nx/2));
        rad=sqrt((RX-Nx/2+.25).^2+(RY-Ny/2+.25).^2);
         u1=-2.5*rad.*(cos(teta*N)+sin(teta*N)).*(Z/Nz)/max(max(max(rad)))+(exp(-((X-Nx/2).^2+(Y-Ny/2).^2+(Z+5).^2)/50));
        %u1=2.5*rad.*(cos(teta*N)+sin(teta*N)).*(RZ/Nz)/max(max(max(rad)))+(exp(-((RX-Nx/2).^2+(RY-Ny/2).^2+(RZ-(R+2)).^2)/20));
%        u1=2.5*(1.5*(exp(-((X-Nx/2).^2+(Y-Ny/2).^2+(Z-(R+2)).^2)/40)));
	%u1 = 2.5*rand(Nx,Ny,Nz);
       %u3=0.-u1;
  	%u1=u1+.2*(rand(Nx,Ny,Nz)-.5);
         u2=(u1+.6*(rand(Nx,Ny,Nz)-.5)).*(0);
         u3=(u1+.6*(rand(Nx,Ny,Nz)-.5)).*(0);
         u4=(u1+.02*(rand(Nx,Ny,Nz)-.5)).*0;
       
% 
%         u(fi<=-.9)=0;
 %        v(fi<=-.9)=0;
%%    store the initial domain

%fi = smooth3(fi,'box',5);
%load temp
fix0(:,:)=fi(Nx/2,:,:);
%%%%%%%%%%% parameters for iteraion loop %%%%%%%%%%%%%%%%%%%%%%%%
step=100;
iter=1;
dt=1e-5;
cont=iter;
ct=0;
%%

    Fm = zeros(Nx,Ny,Nz,NF+1);
    Um = zeros(Nx,Ny,Nz,NF+1);
    Sm = zeros(Nx,Ny,Nz,NF);
    Fm(:,:,:,1) = fi;
    Um(:,:,:,1) = u1;
    
    t = tic();
    disp(1);    
for iter=cont:NF  
    
    for iiter=1:step
 
%%% Phase-field definitons 

        lapfi = lap3Dt(fi);

        gfi=grad3DR(fi);
        gu1=grad3DR(u1);
        gu2=grad3DR(u2);
        gu3=grad3DR(u3);
        gu4=grad3DR(u4);

	gradfi=gfi(:,:,:,1).^2+gfi(:,:,:,2).^2+gfi(:,:,:,3).^2;
	UI=gam1*(gu1(:,:,:,1).*gfi(:,:,:,1)+gu1(:,:,:,2).*gfi(:,:,:,2)+gu1(:,:,:,3).*gfi(:,:,:,3))...
	  +gam2*(gu2(:,:,:,1).*gfi(:,:,:,1)+gu2(:,:,:,2).*gfi(:,:,:,2)+gu2(:,:,:,3).*gfi(:,:,:,3))...
	  +gam3*(gu3(:,:,:,1).*gfi(:,:,:,1)+gu3(:,:,:,2).*gfi(:,:,:,2)+gu3(:,:,:,3).*gfi(:,:,:,3))...
	  +gam4*(gu4(:,:,:,1).*gfi(:,:,:,1)+gu4(:,:,:,2).*gfi(:,:,:,2)+gu4(:,:,:,3).*gfi(:,:,:,3));
	BI=d11*(gu1(:,:,:,1).*gu1(:,:,:,1)+gu1(:,:,:,2).*gu1(:,:,:,2)+gu1(:,:,:,3).*gu1(:,:,:,3))...
	  +d12*(gu1(:,:,:,1).*gu2(:,:,:,1)+gu1(:,:,:,2).*gu2(:,:,:,2)+gu1(:,:,:,3).*gu2(:,:,:,3))...
	  +d13*(gu1(:,:,:,1).*gu3(:,:,:,1)+gu1(:,:,:,2).*gu3(:,:,:,2)+gu1(:,:,:,3).*gu3(:,:,:,3))...
 	  +d14*(gu1(:,:,:,1).*gu4(:,:,:,1)+gu1(:,:,:,2).*gu4(:,:,:,2)+gu1(:,:,:,3).*gu4(:,:,:,3))...
	  +d21*(gu2(:,:,:,1).*gu1(:,:,:,1)+gu2(:,:,:,2).*gu1(:,:,:,2)+gu2(:,:,:,3).*gu1(:,:,:,3))...
	  +d22*(gu2(:,:,:,1).*gu2(:,:,:,1)+gu2(:,:,:,2).*gu2(:,:,:,2)+gu2(:,:,:,3).*gu2(:,:,:,3))...
	  +d23*(gu2(:,:,:,1).*gu3(:,:,:,1)+gu2(:,:,:,2).*gu3(:,:,:,2)+gu2(:,:,:,3).*gu3(:,:,:,3))...
	  +d24*(gu2(:,:,:,1).*gu4(:,:,:,1)+gu2(:,:,:,2).*gu4(:,:,:,2)+gu2(:,:,:,3).*gu4(:,:,:,3))...
	  +d31*(gu3(:,:,:,1).*gu1(:,:,:,1)+gu3(:,:,:,2).*gu1(:,:,:,2)+gu3(:,:,:,3).*gu1(:,:,:,3))...
	  +d32*(gu3(:,:,:,1).*gu2(:,:,:,1)+gu3(:,:,:,2).*gu2(:,:,:,2)+gu3(:,:,:,3).*gu2(:,:,:,3))...
	  +d33*(gu3(:,:,:,1).*gu3(:,:,:,1)+gu3(:,:,:,2).*gu3(:,:,:,2)+gu3(:,:,:,3).*gu3(:,:,:,3))...
	  +d34*(gu3(:,:,:,1).*gu4(:,:,:,1)+gu3(:,:,:,2).*gu4(:,:,:,2)+gu3(:,:,:,3).*gu4(:,:,:,3))...
	  +d41*(gu4(:,:,:,1).*gu1(:,:,:,1)+gu4(:,:,:,2).*gu1(:,:,:,2)+gu4(:,:,:,3).*gu1(:,:,:,3))...
	  +d42*(gu4(:,:,:,1).*gu2(:,:,:,1)+gu4(:,:,:,2).*gu2(:,:,:,2)+gu4(:,:,:,3).*gu2(:,:,:,3))...
	  +d43*(gu4(:,:,:,1).*gu3(:,:,:,1)+gu4(:,:,:,2).*gu3(:,:,:,2)+gu4(:,:,:,3).*gu3(:,:,:,3))...
	  +d44*(gu4(:,:,:,1).*gu4(:,:,:,1)+gu4(:,:,:,2).*gu4(:,:,:,2)+gu4(:,:,:,3).*gu4(:,:,:,3));
	US1=(u1-u11).^2+(u2-u21).^2+(u3-u31).^2+(u4-u41).^2;
	US2=(u1-u12).^2+(u2-u22).^2+(u3-u32).^2+(u4-u42).^2;
	US3=(u1-u13).^2+(u2-u23).^2+(u3-u33).^2+(u4-u43).^2;
	US4=(u1-u14).^2+(u2-u24).^2+(u3-u34).^2+(u4-u44).^2;
	BS=US1.*US2.*US3.*US4;
	BF=(u1-u1f).^2+(u2-u2f).^2+(u3-u3f).^2+(u4-u4f).^2;


%% spontaneous curvature
        Co=bet1*u1.^2+bet2*u2.^2+bet3*u3.^2+bet4*u4.^2;
%% chemistry potential
        mu=((fi-ep1.*Co).*((fi).^2-1)-ep*lapfi);        
%% energy density       
        wfi=0.5*Afi*mu.^2 + 0.5*sigma*gradfi + 0.5*AI*UI + 0.5*As*BI + 0.5*As*(fi.^2 -1).*BS + 0.5*Af*(fi.^2).*BF;


%% Variations

          lapmu = lap3Dt(mu);
%         H=mu;
%         lap3Dt
%         lapmu=lapH;
         
          lapu1 = lap3Dt(u1);
          lapu2 = lap3Dt(u2);
          lapu3 = lap3Dt(u3);
          lapu4 = lap3Dt(u4);

        lapUI=gam1*lapu1+gam2*lapu2+gam3*lapu3+gam4*lapu4;
%% fi
        
        F=Afi*((3*fi.^2-1-2*ep1*fi.*Co).*mu-ep*lapmu)+2*As*fi.*(fi.^2-1).*BS+Af*fi.*BF-0.5*AI*lapUI-sigma*lapfi;

%%  for u's


	Gu1=Afi*(fi.^2-1).*(-2*ep1*bet1*u1)...
	   +As*(fi.^2-1).*((u1-u11).*US2.*US3.*US4+(u1-u12).*US1.*US3.*US4+(u1-u13).*US1.*US2.*US4+(u1-u14).*US1.*US2.*US3)...
	   +Af*(fi.^2).*(u1-u1f)-0.5*As*(2*d11*lapu1+d12*lapu2+d13*lapu3+d14*lapu4)-0.5*AI*gam1*lapfi;
	Gu2=Afi*(fi.^2-1).*(-2*ep1*bet1*u2)...
	   +As*(fi.^2-1).*((u2-u21).*US2.*US3.*US4+(u2-u22).*US1.*US3.*US4+(u2-u23).*US1.*US2.*US4+(u2-u24).*US1.*US2.*US3)...
	   +Af*(fi.^2).*(u2-u2f)-0.5*As*(d21*lapu1+2*d22*lapu2+d23*lapu3+d24*lapu4)-0.5*AI*gam2*lapfi;
	Gu3=Afi*(fi.^2-1).*(-2*ep1*bet1*u3)...
	   +As*(fi.^2-1).*((u3-u31).*US2.*US3.*US4+(u3-u32).*US1.*US3.*US4+(u3-u33).*US1.*US2.*US4+(u3-u34).*US1.*US2.*US3)...
	   +Af*(fi.^2).*(u3-u3f)-0.5*As*(d31*lapu1+d32*lapu2+2*d33*lapu3+d34*lapu4)-0.5*AI*gam3*lapfi;
	Gu4=Afi*(fi.^2-1).*(-2*ep1*bet1*u4)...
	   +As*(fi.^2-1).*((u4-u41).*US2.*US3.*US4+(u4-u42).*US1.*US3.*US4+(u4-u43).*US1.*US2.*US4+(u4-u44).*US1.*US2.*US3)...
	   +Af*(fi.^2).*(u4-u4f)-0.5*As*(d41*lapu1+d42*lapu2+d43*lapu3+2*d44*lapu4)-0.5*AI*gam4*lapfi;
%%%


	gmu=grad3DR(mu);
	ggfi1=grad3DR(gfi(:,:,:,1));
	ggfi2=grad3DR(gfi(:,:,:,2));        
	ggfi3=grad3DR(gfi(:,:,:,3));  

%% Stress tensor


	str(:,:,:,1,1)=(wfi(:,:,:) -fi(:,:,:).*F(:,:,:)-sigma*gfi(:,:,:,1).*gfi(:,:,:,1)...
	-0.5*AI*(gam1*gu1(:,:,:,1)+gam2*gu2(:,:,:,1)+gam3*gu3(:,:,:,1)+gam4*gu4(:,:,:,1)).*gfi(:,:,:,1)-Afi*ep*gmu(:,:,:,1).*gfi(:,:,:,1)+Afi*ep*mu.*ggfi1(:,:,:,1));
    str(:,:,:,2,2)=(wfi(:,:,:) -fi(:,:,:).*F(:,:,:)-sigma*gfi(:,:,:,2).*gfi(:,:,:,2)...
	-0.5*AI*(gam1*gu1(:,:,:,2)+gam2*gu2(:,:,:,2)+gam3*gu3(:,:,:,2)+gam4*gu4(:,:,:,2)).*gfi(:,:,:,2)-Afi*ep*gmu(:,:,:,2).*gfi(:,:,:,2)+Afi*ep*mu.*ggfi2(:,:,:,2));
    str(:,:,:,3,3)=(wfi(:,:,:) -fi(:,:,:).*F(:,:,:)-sigma*gfi(:,:,:,3).*gfi(:,:,:,3)...
	-0.5*AI*(gam1*gu1(:,:,:,3)+gam2*gu2(:,:,:,3)+gam3*gu3(:,:,:,3)+gam4*gu4(:,:,:,3)).*gfi(:,:,:,3)-Afi*ep*gmu(:,:,:,3).*gfi(:,:,:,3)+Afi*ep*mu.*ggfi3(:,:,:,3));
 
	str(:,:,:,1,2)=(-sigma*gfi(:,:,:,1).*gfi(:,:,:,2)...
	-0.5*AI*(gam1*gu1(:,:,:,2)+gam2*gu2(:,:,:,2)+gam3*gu3(:,:,:,2)+gam4*gu4(:,:,:,2)).*gfi(:,:,:,1)-Afi*ep*gmu(:,:,:,2).*gfi(:,:,:,1)+Afi*ep*mu.*ggfi2(:,:,:,1));
	str(:,:,:,1,3)=(-sigma*gfi(:,:,:,1).*gfi(:,:,:,3)...
	-0.5*AI*(gam1*gu1(:,:,:,3)+gam2*gu2(:,:,:,3)+gam3*gu3(:,:,:,3)+gam4*gu4(:,:,:,3)).*gfi(:,:,:,1)-Afi*ep*gmu(:,:,:,3).*gfi(:,:,:,1)+Afi*ep*mu.*ggfi3(:,:,:,1));
	str(:,:,:,2,1)=(-sigma*gfi(:,:,:,2).*gfi(:,:,:,1)...
	-0.5*AI*(gam1*gu1(:,:,:,1)+gam2*gu2(:,:,:,1)+gam3*gu3(:,:,:,1)+gam4*gu4(:,:,:,1)).*gfi(:,:,:,2)-Afi*ep*gmu(:,:,:,1).*gfi(:,:,:,2)+Afi*ep*mu.*ggfi1(:,:,:,2));
	str(:,:,:,2,3)=(-sigma*gfi(:,:,:,2).*gfi(:,:,:,3)...
	-0.5*AI*(gam1*gu1(:,:,:,3)+gam2*gu2(:,:,:,3)+gam3*gu3(:,:,:,3)+gam4*gu4(:,:,:,3)).*gfi(:,:,:,2)-Afi*ep*gmu(:,:,:,3).*gfi(:,:,:,2)+Afi*ep*mu.*ggfi3(:,:,:,2));
	str(:,:,:,3,1)=(-sigma*gfi(:,:,:,3).*gfi(:,:,:,1)...
	-0.5*AI*(gam1*gu1(:,:,:,1)+gam2*gu2(:,:,:,1)+gam3*gu3(:,:,:,1)+gam4*gu4(:,:,:,1)).*gfi(:,:,:,3)-Afi*ep*gmu(:,:,:,1).*gfi(:,:,:,3)+Afi*ep*mu.*ggfi1(:,:,:,3));
	str(:,:,:,3,2)=(-sigma*gfi(:,:,:,3).*gfi(:,:,:,2)...
	-0.5*AI*(gam1*gu1(:,:,:,2)+gam2*gu2(:,:,:,2)+gam3*gu3(:,:,:,2)+gam4*gu4(:,:,:,2)).*gfi(:,:,:,3)-Afi*ep*gmu(:,:,:,2).*gfi(:,:,:,3)+Afi*ep*mu.*ggfi2(:,:,:,3)); 


%% 
%% phase-field dynamics

  	gGu1=grad3DR(Gu1);
	gGu2=grad3DR(Gu2);
	gGu3=grad3DR(Gu3);
	gGu4=grad3DR(Gu4);
	% this is the proyection of stress tensor and the variations of u
 	du1S1(:,:,:)=str(:,:,:,1,1).*gGu1(:,:,:,1)+str(:,:,:,1,2).*gGu1(:,:,:,2)+str(:,:,:,1,3).*gGu1(:,:,:,3); 
 	du1S2(:,:,:)=str(:,:,:,2,1).*gGu1(:,:,:,1)+str(:,:,:,2,2).*gGu1(:,:,:,2)+str(:,:,:,2,3).*gGu1(:,:,:,3);
 	du1S3(:,:,:)=str(:,:,:,3,1).*gGu1(:,:,:,1)+str(:,:,:,3,2).*gGu1(:,:,:,2)+str(:,:,:,3,3).*gGu1(:,:,:,3);
 	du2S1(:,:,:)=str(:,:,:,1,1).*gGu2(:,:,:,1)+str(:,:,:,1,2).*gGu2(:,:,:,2)+str(:,:,:,1,3).*gGu2(:,:,:,3); 
 	du2S2(:,:,:)=str(:,:,:,2,1).*gGu2(:,:,:,1)+str(:,:,:,2,2).*gGu2(:,:,:,2)+str(:,:,:,2,3).*gGu2(:,:,:,3);
 	du2S3(:,:,:)=str(:,:,:,3,1).*gGu2(:,:,:,1)+str(:,:,:,3,2).*gGu2(:,:,:,2)+str(:,:,:,3,3).*gGu2(:,:,:,3);
 	du3S1(:,:,:)=str(:,:,:,1,1).*gGu3(:,:,:,1)+str(:,:,:,1,2).*gGu3(:,:,:,2)+str(:,:,:,1,3).*gGu3(:,:,:,3); 
 	du3S2(:,:,:)=str(:,:,:,2,1).*gGu3(:,:,:,1)+str(:,:,:,2,2).*gGu3(:,:,:,2)+str(:,:,:,2,3).*gGu3(:,:,:,3);
 	du3S3(:,:,:)=str(:,:,:,3,1).*gGu3(:,:,:,1)+str(:,:,:,3,2).*gGu3(:,:,:,2)+str(:,:,:,3,3).*gGu3(:,:,:,3);
 	du4S1(:,:,:)=str(:,:,:,1,1).*gGu4(:,:,:,1)+str(:,:,:,1,2).*gGu4(:,:,:,2)+str(:,:,:,1,3).*gGu4(:,:,:,3); 
 	du4S2(:,:,:)=str(:,:,:,2,1).*gGu4(:,:,:,1)+str(:,:,:,2,2).*gGu4(:,:,:,2)+str(:,:,:,2,3).*gGu4(:,:,:,3);
 	du4S3(:,:,:)=str(:,:,:,3,1).*gGu4(:,:,:,1)+str(:,:,:,3,2).*gGu4(:,:,:,2)+str(:,:,:,3,3).*gGu4(:,:,:,3);

	 %this is the divergence
 	gu1s1=grad3DR(du1S1(:,:,:));
 	gu1s2=grad3DR(du1S2(:,:,:));
 	gu1s3=grad3DR(du1S3(:,:,:));
 	gu2s1=grad3DR(du2S1(:,:,:));
 	gu2s2=grad3DR(du2S2(:,:,:));
 	gu2s3=grad3DR(du2S3(:,:,:));
 	gu3s1=grad3DR(du3S1(:,:,:));
 	gu3s2=grad3DR(du3S2(:,:,:));
 	gu3s3=grad3DR(du3S3(:,:,:));
 	gu4s1=grad3DR(du4S1(:,:,:));
 	gu4s2=grad3DR(du4S2(:,:,:));
 	gu4s3=grad3DR(du4S3(:,:,:));


 	Su1=gu1s1(:,:,:,1)+gu1s2(:,:,:,2)+gu1s3(:,:,:,3);
 	Su2=gu2s1(:,:,:,1)+gu2s2(:,:,:,2)+gu2s3(:,:,:,3);
 	Su3=gu3s1(:,:,:,1)+gu3s2(:,:,:,2)+gu3s3(:,:,:,3);
 	Su4=gu4s1(:,:,:,1)+gu4s2(:,:,:,2)+gu4s3(:,:,:,3);
 
%%            dynamical equations,  for conservation of mass of fi and u's     

          lapF = lap3Dt(F);       
      
          lapGu1 = lap3Dt(Gu1);
          lapGu2 = lap3Dt(Gu2);
          lapGu3 = lap3Dt(Gu3);
          lapGu4 = lap3Dt(Gu4);          

         %I=200.*(u1+u2+u3+u4);
         I=alpha.*Gu1;
         I(find(abs(fi)>=.9))=0;
         
        fi=fi+Dfi*dt*(lapF+I);
        fi(:,:,1)=fi(:,:,2);
       % Iu=sum(sum(sum(u)))/Nx/Ny/Nz-u0;
        u1=u1+dt*(Du1*lapGu1+eta*Su1);
        u2=u2+dt*(Du2*lapGu2+eta*Su2);
        u3=u3+dt*(Du3*lapGu3+eta*Su3);
        u4=u4+dt*(Du4*lapGu4+eta*Su4);

        %uts
       % u(:,:,1)=u(:,:,2);  
        noFlux2(fi,u1); 
    end
    u1x(:,:)=u1(:,Ny/2,:);
    fix(:,:)=fi(Nx/2,:,:);

    %%   
hh=max(max(max(isnan(fi(:,:,:)))));
    if hh==1;
       disp(NaN);
        break
    end

%para hacer peliculas 3D usar cine3D y cine2D

    Fm(:,:,:,iter+1)=fi(:,:,:);
    Um(:,:,:,iter+1)=u1(:,:,:);
    Sm(:,:,:,iter)=Su1(:,:,:);


    disp(iter+1)%time loop
end

time = toc(t);

save('julio1c.mat');
