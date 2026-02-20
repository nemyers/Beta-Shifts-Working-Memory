xyz = randn(3,12);

figure
hold on

%plot dot cloud
plot3(xyz(1,:),xyz(2,:),xyz(3,:),'ko')
%plot projections
projection_constant = -5;
%projection onto X
plot3(xyz(1,:)*0+projection_constant,xyz(2,:),xyz(3,:),'bo')
%projection onto Y
plot3(xyz(1,:),xyz(2,:)*0+projection_constant,xyz(3,:),'go')
%projection onto Z
plot3(xyz(1,:),xyz(2,:),xyz(3,:)*0 + projection_constant,'ro')

set(gca,'cameraposition',[1 -1 1],'YGrid','on','XGrid','on','ZGrid','on')