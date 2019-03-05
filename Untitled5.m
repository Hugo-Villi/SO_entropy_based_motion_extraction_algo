function 3D_scatter_plot = scatter_3D_plot(markers_values)
k=1;
for i=1:3:size(markers_values,2)
    x(:,k)=markers_values(:,i);
    y(:,k)=markers_values(:,i+1);
    z(:,k)=markers_values(:,i+2);
    k=k+1;
end
for i=1:size(x,1)
scatter3(x(50,i),y(50,i),z(50,i))


end