%%acquisition of the needed data from the motion capture.
file = uigetfile('*.c3d');
acq = btkReadAcquisition(file); %read the file
markers = btkGetMarkers(acq);   %get the information of the markers in a structure
labels = fieldnames(markers);   %get the labels of the markers
markers_values = btkGetMarkersValues(acq); %give an array filled with the coordinates of the marker (x,y,z of first marker, x,y,z of second, and so on)

%%calculation of the displacement vector
num_row = size(markers_values,1);   %a row = a frame
num_col = size(markers_values,2);   %this is equal to 3*number of markers
displacement=zeros(num_row,num_col);    %preallocation
for i=1:(num_row-1) %-1 otherwise it will crash trying to reach the last frame + 1
    for j=1:num_col
        displacement(i,j)=markers_values(i+1,j)-markers_values(i,j);    %simple soustraction between frame+1-frame
    end
end
%give only displacements for x, y and z, to ease the next steps (max/min
%detection and discretization)
for i=3:3:num_col
    displacement_x(:,i/3)=displacement(:,i-2);
end
for i=3:3:num_col
    displacement_y(:,i/3)=displacement(:,i-1);
end
for i=3:3:num_col
    displacement_z(:,i/3)=displacement(:,i);
end
%%displacement histogram
n = 12; %number of discretizations levels, this setting may have importance on the results
maximum=max(displacement);  %gives the max values for all the column hence for all the x,y,z coordinates respectively
minimum=min(displacement);  %same for min
%the following part get the maximum and minimums of x,y and z displacement
%to set the range of the histograms
max_x = max(displacement_x,[],'all');
min_x = min(displacement_x,[],'all');

max_y = max(displacement_y,[],'all');
min_y = min(displacement_y,[],'all');

max_z = max(displacement_z,[],'all');
min_z = min(displacement_z,[],'all');

%creating the histograms
discretization_x=linspace(min_x,max_x,n);   %the function linspace creates a vector of values envenly distributed along a range
for i=1:size(displacement_x,1)    %will generate an histogram for each frame for the x coordinates
    histogram_x(i,:)=histcounts(displacement_x(i,:),discretization_x);
end
discretization_y=linspace(min_y,max_y,n);
for i=1:size(displacement_y,1)
    histogram_y(i,:)=histcounts(displacement_y(i,:),discretization_y);
end
discretization_z=linspace(min_z,max_z,n);
for i=1:size(displacement_z,1)
    histogram_z(i,:)=histcounts(displacement_z(i,:),discretization_z);
end


%computing the probabilities
for i=1:size(histogram_x)
    for j=1:size(histogram_x,2)
        prob_x(i,j)=histogram_x(i,j)/(num_col/3);
    end
end
for i=1:size(histogram_x)
    for j=1:size(histogram_x,2)
        prob_y(i,j)=histogram_y(i,j)/(num_col/3);
    end
end
for i=1:size(histogram_x)
    for j=1:size(histogram_x,2)
        prob_z(i,j)=histogram_z(i,j)/(num_col/3);
    end
end


%filling the C(r,s) matrix
%with the same global range as the histograms
%create the matrices for each marker
[matrix_r_s_global_x,matrix_r_s_x,Cx,proba_r_s_x]=compute_C(histogram_x,displacement_x,discretization_x,labels);
[matrix_r_s_global_y,matrix_r_s_y,Cy,proba_r_s_y]=compute_C(histogram_y,displacement_y,discretization_y,labels);
[matrix_r_s_global_z,matrix_r_s_z,Cz,proba_r_s_z]=compute_C(histogram_z,displacement_z,discretization_z,labels);

Ix=mutual_info(displacement,Cx,prob_x);
Iy=mutual_info(displacement,Cy,prob_y);
Iz=mutual_info(displacement,Cz,prob_z);
I=Ix+Iy+Iz;


