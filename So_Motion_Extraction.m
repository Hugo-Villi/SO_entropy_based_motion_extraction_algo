%%acquisition of the needed data from the motion capture.
clearvars;
file = uigetfile('*.c3d');                  %open a window for the user to get the file
acq = btkReadAcquisition(file);             %read the file
labels = fieldnames(btkGetMarkers(acq));    %get the labels of the markers
markers_values = btkGetMarkersValues(acq);  %give an array filled with the coordinates of the marker (x,y,z of first marker, x,y,z of second, and so on)

%%calculation of the displacement vector
displacement=zeros(size(markers_values,1)-1,size(markers_values,2));%preallocation
for i=1:(size(markers_values,1)-1)                                  %-1 otherwise it will crash trying to reach the last frame + 1
    displacement(i,:)=markers_values(i+1,:)-markers_values(i,:);    %simple soustraction between frame+1-frame
end
%give only displacements for x, y and z, to ease the next steps (max/min
%detection and discretization)
for i=3:3:size(markers_values,2)
    displacement_x(:,i/3)=displacement(:,i-2);
    displacement_y(:,i/3)=displacement(:,i-1);
    displacement_z(:,i/3)=displacement(:,i);
end
%%displacement histogram
maximum=max(displacement,[],'all');     %gives the max values for all the column hence for all the x,y,z coordinates respectively
minimum=min(displacement,[],'all');     %same for min
%the following part get the maximum and minimums of x,y and z displacement
%to set the range of the histograms
max_x = max(displacement_x,[],'all');
min_x = min(displacement_x,[],'all');

max_y = max(displacement_y,[],'all');
min_y = min(displacement_y,[],'all');

max_z = max(displacement_z,[],'all');
min_z = min(displacement_z,[],'all');

%creating the histograms
n = 9;                                %number of discretizations levels, this setting may have importance on the results
discretization_x=linspace(min_x,max_x,n);   %the function linspace creates a vector of values envenly distributed along a range
discretization_y=linspace(min_y,max_y,n);
discretization_z=linspace(min_z,max_z,n);
discretization_global=linspace(minimum,maximum,n);
for i=1:size(displacement_x,1)    %will generate an histogram for each frame for the x,y and z coordinates
    histogram_x(i,:)=histcounts(displacement_x(i,:),discretization_x);
    histogram_y(i,:)=histcounts(displacement_y(i,:),discretization_y);
    histogram_z(i,:)=histcounts(displacement_z(i,:),discretization_z);
    histogram_max_min(i,:)=histcounts(displacement(i,:),discretization_global);
end

%computing the probabilities, simply bu dividing the histogram by the
%number of marker
for i=1:size(histogram_x)
        prob_x(i,:)=histogram_x(i,:)/(size(markers_values,2)/3);
        prob_y(i,:)=histogram_y(i,:)/(size(markers_values,2)/3);
        prob_z(i,:)=histogram_z(i,:)/(size(markers_values,2)/3);
        prob_global(i,:)=histogram_max_min(i,:)/(size(markers_values,2)/3);    
end

%computation of the C matrix
Cx=compute_C(n,displacement_x,discretization_x,labels);
Cy=compute_C(n,displacement_y,discretization_y,labels);
Cz=compute_C(n,displacement_z,discretization_z,labels);
C_max_min=compute_C(n,displacement,discretization_global,labels);
%compute the mutual information for each dimension and sum it to get the
%total mutual information I
Ix=mutual_info(size(displacement,1),Cx,prob_x);
Iy=mutual_info(size(displacement,1),Cy,prob_y);
Iz=mutual_info(size(displacement,1),Cz,prob_z);
I=Ix+Iy+Iz;
I_max_min=mutual_info(size(displacement,1),C_max_min,prob_global);

plot(I,'marker','o');
hold on
plot(Ix,'marker','+');
hold on
plot(Iy,'marker','*');
hold on
plot(Iz,'marker','x');
hold on
plot(I_max_min,'marker','s');
hold on
xlabel('frames') 
ylabel('Mutual information I') 
legend('I','Ix','Iy','Iz','I_max_min');
hold off
%settings for the keypose detection
window_size=60;
threshold=-1;
%execute the function to detect the keyposes
[keyposes_temp,I_localized,I_trimmed]=keyposes_detection(I,window_size,threshold);