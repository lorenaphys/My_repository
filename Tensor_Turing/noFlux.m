function noFlux(f,g)

h = size(f);
Nx = h(1);
Ny = h(2);
Nz = h(3);

% for i = 1:Nx
%    for j = 1:Ny
%       for k = 1:Nz
%          if abs(f(i,j,k)) <= 0.1
%              if i == 1
%                  g(i+1,j,k) = g(i,j,k);
%              else
%                  g(i-1,j,k) = g(i,j,k);
%              end
%              
%              if j == 1
%                  g(i,j+1,k) = g(i,j,k);
%              else
%                  g(i,j-1,k) = g(i,j,k);
%              end
%              
%              if k == 1
%                  g(i,j,k+1) = g(i,j,k);
%              else
%                  g(i,j,k-1) = g(i,j,k);
%              end
%          end
%       end
%    end
% end

for i = 1:Nx
   for j = 1:Ny
      for k = 1:Nz
          if abs(f(i,j,k)) <= 0.1
            if i < Nx
			a = f(i+1,j,k);
		    else
			a = f(i-1,j,k);
            end
            b = g(i,j,k);
            if j < Ny
                c = f(i,j+1,k);
                d = g(i,j+1,k);
            else
                c = f(i,j-1,k);
                d = g(i,j-1,k);
            end
            if k < Nz
                e = f(i,j,k+1);
                h = g(i,j,k+1);
            else
                e = f(i,j,k-1);
                h = g(i,j,k-1); 
            end
            g(i+1,j,k) = (b-c*(d-b)-e*(h-b))/a;
          end
      end
   end
end
