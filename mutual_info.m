function Ix = mutual_info(frame_number,Cx,prob_x)
for k=1:frame_number  %go through each frame
    sum=0;  
    for i=1:size(Cx,1)  %go through the row and column of the n*n C(r,s) matrices
        for j=1:size(Cx,2)
            if Cx(i,j)~=0 && prob_x(k,i)~=0 && prob_x(k,j)  %if C or prob is equal to 0 the resulting value would be infinite and the zeros will therefore be excluded
                sum=sum+(-Cx(i,j)*log(Cx(i,j)/(prob_x(k,i)*prob_x(k,j))));  %the equation given in So (2003) to get the mutual information measure
            end
        end
    end
    Ix(k)=sum;      
end