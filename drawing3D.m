function drawing3D(R1,a)
[x,y,z] = meshgrid(-10*a:10*a/100:10*a);
c1 =-R1^2+(y.^2+(x-R1+a).^2+z.^2);
c2 = -R1^2+y.^2+(x+R1-a).^2+z.^2;
c3=zeros(size(c1));
c = max(max(c1,c2),c3);
c11 = c1<=0 & c2<=0 & c3<=0;
fv1 = isosurface(x,y,z,c,0);
fv2 = isosurface(x,y,z,c11,0);
p = patch(fv1);
isonormals(x,y,z,c,p)
set(p,'facecolor',[0 .5 1],'edgecolor','none')
view(150,30),axis image,grid on
camlight
lighting gouraud
figure
p1 = patch(fv2);
isonormals(x,y,z,c11,p1)
set(p1,'facecolor',[0 .5 1],'edgecolor','none')
view(150,30),axis image,grid on
camlight
lighting gouraud