function noFlux(f)

h = size(f);
Nx = h(1);
Ny = h(2);
Nz = h(3);

for i = 2:Nx
   for j = 2:Ny
      for k = 2:Nz
         if abs(f(i,j,k)) <= 0.1
             f(i,j,k) = f(i-1,j-1,k-1); 
         end
      end
   end
end
