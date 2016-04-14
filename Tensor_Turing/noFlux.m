function noFlux(f)

h = size(f);
Nx = h(1);
Ny = h(2);
Nz = h(3);

for i = 1:Nx
   for j = 1:Ny
      for k = 1:Nz
        if f == 0
           if (i < Nx) || (j < Ny) || (k < Nz) 
               f(i,j,k) = f(i+1,j+1,k+1);
           elseif i == Nx
               f(i,j,k) = f(i-1,j,k);
           elseif j == Ny
               f(i,j,k) = f(i,j-1,k);
           else
               f(i,j,k) = f(i,j,k-1);
           end
        end
      end
   end
end
