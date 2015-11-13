load('oct21b.mat')
f = 41;
fi(:,:,:) = Fm(:,:,:,f);
v = Vm(:,1);

fi = smooth3(fi,'box',5);
isosurface(fi,0);
axis equal, view(-14,40), axis off
colormap jet 
title(['Volumen reducido, v = ',num2str(v(f))])
