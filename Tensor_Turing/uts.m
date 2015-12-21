% Program to calculate turing BVAM 3 dimensions

 
 
 %   Laplacoano con grad=0 en la frontera del cubo 

     
 %for ii=1:5
         H=u;
        lap3Dt
        lapu=lapH;
        

        H=v;
        lap3Dt
        lapv=lapH; 
%         
        lapu(fi<=-.9)=0;
        lapv(fi<=-.9)=0;
        u(fi<=-.9)=0;
        v(fi<=-.9)=0;

%       u=u+dt*(Du*lapu+(u+a*v-c*u.*v-u.*v.^2));

        u=u+dt1*(Du*(lapu)+(u+a*v-c*u.*v-u.*v.^2));
        v=v+dt1*(Dv*lapv+(b*v+h*u+c*u.*v+u.*v.^2));        
%        u(:,:,1)=u(:,:,2);
%       v(:,:,1)=v(:,:,2);
 %end
 

 
