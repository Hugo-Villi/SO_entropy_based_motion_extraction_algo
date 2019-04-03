function proba_r_s_global=compute_C(n,displacement,discretization,labels)
total_number_of_event=(size(displacement,1)-1)*size(labels,1);    %get the total number of events to get the global probabilities
matrix_r_s_global=zeros(n-1,n-1);   %declare the matrix

for i=1:size(labels,1)
    for j=1:size(displacement,1)-1
        r=find(histcounts(displacement(j,i),discretization)==1);    %get the position in the histogram at the frame f
        s=find(histcounts(displacement(j+1,i),discretization)==1);  %get the position in the histogram at the frame f+1
        matrix_r_s_global(r,s)=matrix_r_s_global(r,s)+1;    %same but for the global matrix        
    end
end

proba_r_s_global=matrix_r_s_global/total_number_of_event;   %compute the global proba
end