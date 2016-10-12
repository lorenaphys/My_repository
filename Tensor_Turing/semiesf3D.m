Nx=50;
Ny=50;
Nz=70;
R=15;

% load initial
% u0=uout+ues; v0=vout+ves; w0=wout+wes; x0=xout+xes;
fi=ones(Nx,Ny,Nz);
r = zeros(Nx,Ny,Nz);

ancho=2;
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
      r(i,j,k)=sqrt((i-Nx/2)^2+(j-Ny/2)^2+(k)^2);
      if r(i,j,k)>=R
      fi(i,j,k)=-1;
      end
        end
   
   end
end

fi = smooth3(fi,'box',5);
isosurface(fi,0)
axis equal, view (-16,24), axis off
axis([1 50 1 50 1 70]),
