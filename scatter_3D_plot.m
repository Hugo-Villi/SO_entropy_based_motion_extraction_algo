function scatter_plot = scatter_3D_plot(markers_values,j)
k=1;
for i=1:3:size(markers_values,2)
    x(:,k)=markers_values(:,i);
    y(:,k)=markers_values(:,i+1);
    z(:,k)=markers_values(:,i+2);
    k=k+1;
end

for i=1:size(x,2)
    scatter_plot=scatter3(x(j,i),y(j,i),z(j,i));
    axis('equal');
    %view(2)
    hold on
end

hold off
end