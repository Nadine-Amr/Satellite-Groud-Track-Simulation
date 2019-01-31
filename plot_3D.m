function [] = plot_3D(r,v)

for i=1:length(r)
    x(i) = r{i}(1);
    y(i) = r{i}(2);
    z(i) = r{i}(3);
end
plot3(x, y, z)
hold on
[xx, yy, zz] = sphere(100);
surf(6378*xx, 6378*yy, 6378*zz)