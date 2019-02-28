function [matrix_r_s_global,matrix_r_s,proba_r_s_global,proba_r_s]=compute_C(histogram_x,displacement_x,discretization_x,labels)
total_number_of_event=0;    %get the total number of events to get the global probabilities
matrix_r_s_global=zeros(size(histogram_x,2),size(histogram_x,2));   %declare the matrix
for i=1:size(labels,1)  %create the structure to compute the proba for each marker
    matrix_r_s.(labels{i})=zeros(size(histogram_x,2),size(histogram_x,2));
    proba_r_s.(labels{i})=zeros(size(histogram_x,2),size(histogram_x,2));
end
for i=1:size(displacement_x,2)
    for j=1:size(displacement_x,1)-1
        r=find(histcounts(displacement_x(j,i),discretization_x)==1);    %get the position in the histogram at the frame f
        s=find(histcounts(displacement_x(j+1,i),discretization_x)==1);  %get the position in the histogram at the frame f+1
        matrix_r_s.(labels{i})(r,s)=matrix_r_s.(labels{i})(r,s)+1;  %add a 1 in the matrix to store the change
        matrix_r_s_global(r,s)=matrix_r_s_global(r,s)+1;    %same but for the global matrix
        total_number_of_event=total_number_of_event+1;  
    end
end
for i=1:size(labels,1)  %compute the proba for each markers
    proba_r_s.(labels{i})=matrix_r_s.(labels{i})/size(displacement_x,1);
end
proba_r_s_global=matrix_r_s_global/total_number_of_event;   %compute the global proba
end