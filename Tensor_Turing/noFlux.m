function noFlux(f)

h = size(f);
Nx = h(1);
Ny = h(2);
Nz = h(3);

for i = 1:Nx
   for j = 1:Ny
      for k = 1:Nz
         if abs(f(i,j,k)) <= 0.1
             if i == 1
                 f(i+1,j,k) = f(i,j,k);
             else
                 f(i-1,j,k) = f(i,j,k);
             end
             
             if j == 1
                 f(i,j+1,k) = f(i,j,k);
             else
                 f(i,j-1,k) = f(i,j,k);
             end
             
             if k == 1
                 f(i,j,k+1) = f(i,j,k);
             else
                 f(i,j,k-1) = f(i,j,k);
             end
         end
      end
   end
end
