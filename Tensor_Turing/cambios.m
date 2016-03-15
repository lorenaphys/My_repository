function cambios(X,X1)

x = size(X);
y = size(x);
a = 0;

if y(2) == 3
Nx = x(1);
Ny = x(2);
Nz = x(3);

e = zeros(Nx,Ny,Nz);

for i = 1:Nx
    for j = 1:Ny
       for k = 1:Nz
          e(i,j,k) = abs(X(i,j,k)-X1(i,j,k));
          if X(i,j,k)~= X1(i,j,k)
             a = 1;
          end
       end
    end
end

p = mean(mean(mean(e)));

end

if y == 2
   Nx = x(1);
   Ny = x(2);
   e = zeros(Nx,Ny);
   
   for i = 1:Nx
      for j = 1:Ny
          e(i,j) = abs(X(i,j)-X1(i,j));
          if X(i,j) ~= X1(i,j)
             a = 1; 
          end
      end
   end
 p = mean(mean(e));  
   
end


if a == 0
    disp('No hubo cambios')
else
    fprintf('Hubo cambios, con una desviacion en promedio de %d \n',p)
end