function[gradH]=grad3DR(H)

dimen=size(H);

Nx=dimen(1);
Ny=dimen(2);
Nz=dimen(3);

    
grad=0*H;
gradUx=grad;
gradUy=grad;
gradUz=grad;

gradUx(2:Nx-1,:,:)=(H(3:Nx,:,:)-H(1:Nx-2,:,:))/2;
gradUy(:,2:Ny-1,:)=(H(:,3:Ny,:)-H(:,1:Ny-2,:))/2;
gradUz(:,:,2:Nz-1)=(H(:,:,3:Nz)-H(:,:,1:Nz-2))/2;

     gradUy(:,1,:)=0;
     gradUy(:,2,:)=0;
     gradUy(:,Ny,:)=0;
     gradUy(:,Ny-1,:)=0;
     gradUx(1,:,:)=0;
     gradUx(2,:,:)=0;     
     gradUx(Nx,:,:)=0;
     gradUx(Nx-1,:,:)=0;
     gradUz(:,:,1)=0;
     gradUz(:,:,2)=0;  
     gradUz(:,:,Nz)=0;
     gradUz(:,:,Nz-1)=0;
          

 gradH=zeros(Nx,Ny,Nz,3);    
     
     gradH(:,:,:,1)=gradUx;
     gradH(:,:,:,2)=gradUy;
     gradH(:,:,:,3)=gradUz;
     
