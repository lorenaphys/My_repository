aux=H;

    

lap=0*aux;
LapUx=lap;
LapUy=lap;
LapUz=lap;


for n=2:Ny-1
      LapUy(:,n,:)=aux(:,n-1,:)-2*aux(:,n,:)+aux(:,n+1,:); 
end
 for m=2:Nx-1;
      LapUx(m,:,:)=aux(m-1,:,:)-2*aux(m,:,:)+aux(m+1,:,:);
 end
for mm=2:Nz-1;
      LapUz(:,:,mm)=aux(:,:,mm-1)-2*aux(:,:,mm)+aux(:,:,mm+1);
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