function[lapH]=lapf3D(fun)

dimen=size(fun);

Nx=dimen(1);
Ny=dimen(2);
Nz=dimen(3);
%fun=H;

    

lap=gpuArray(zeros(Nx,Ny,Nz));
LapUx=lap;
LapUy=lap;
LapUz=lap;


for n=2:Ny-1
      LapUy(:,n,:)=fun(:,n-1,:)-2*fun(:,n,:)+fun(:,n+1,:); 
end
 for m=2:Nx-1;
      LapUx(m,:,:)=fun(m-1,:,:)-2*fun(m,:,:)+fun(m+1,:,:);
 end
for mm=2:Nz-1;
      LapUz(:,:,mm)=fun(:,:,mm-1)-2*fun(:,:,mm)+fun(:,:,mm+1);
end
     LapUy(:,1,:)=0;%aux(2,:)-aux(1,:);
     LapUy(:,Ny,:)=0;%aux(Nx-1,:)-aux(Nx,:);
     LapUx(1,:,:)=0;%aux(:,2)-aux(:,1);
     LapUx(Nx,:,:)=0;%aux(:,Ny-1)-aux(:,Ny);
     LapUz(:,:,1)=0;
     LapUz(:,:,Nz)=0;
     lapH=LapUx+LapUy+LapUz;
     
% for qq=1:Nx
%     for rr=1:Ny
%         for ss=1:Nz
%         lapH(qq,rr,ss)=lap(qq,rr,ss);
%         end
%     end
% end