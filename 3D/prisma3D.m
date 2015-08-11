
fi=-1*ones(Nx,Ny,Nz);

ancho=2;
for i=20:60
    for j=20:30
        for k=20:30
      
        fi(i,j,k)=1;
     
        end
   
   end
end

%fi = smooth3(fi,'box',5);

%isosurface(fi,0)
