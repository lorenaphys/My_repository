function revision(f)

h = size(f);
Nx = h(1);
Ny = h(2);
Nz = h(3);

for i = 1:Nx
   for j = 1:Ny
      for k = 1:Nz
         a = isnan(f(i,j,k));
         if a == 1
            break
         end
      end
   end
end