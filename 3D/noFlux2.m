function noFlux2(f,g)

h = size(f);
Nx = h(1);
Ny = h(2);
Nz = h(3);

for i = 1:Nx
   for j = 1:Ny
      for k = 1:Nz
         if abs(f(i,j,k))<0.1
            b = g(i-1,j,k);
            if j == 1
                c = (f(i,j+1,k)-f(i,j,k))*(g(i,j+1,k)-g(i,j,k));
            elseif j > 1 && j < Ny
                c = (f(i,j+1,k)-f(i,j-1,k))*(g(i,j+1,k)-g(i,j-1,k));
            else
                c = (f(i,j,k)-f(i,j-1,k))*(g(i,j,k)-g(i,j-1,k));
            end
            if k == 1
               d = (f(i,j,k+1)-f(i,j,k))*(g(i,j,k+1)-g(i,j,k));
            elseif k > 1 && k < Nz
               d = (f(i,j,k+1)-f(i,j,k-1))*(g(i,j,k+1)-g(i,j,k-1));
            else
               d = (f(i,j,k)-f(i,j,k-1))*(g(i,j,k)-g(i,j,k-1));
            end
            if i == 1
               a = f(i+1,j,k)-f(i,j,k);
               g(i+1,j,k) = (b+c+d)/a;
            elseif i > 1 && i < Nx
               a = f(i+1,j,k)-f(i-1,j,k);
               g(i+1,j,k) = (b+c+d)/a;
            else
               a = f(i,j,k)-f(i-1,j,k);
               g(i-1,j,k) = -(b+c+d)/a;
            end
         end
      end
   end
end