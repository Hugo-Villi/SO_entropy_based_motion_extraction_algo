function Ix = mutual_info(displacement,Cx,prob_x)
for k=1:size(displacement,1)
    sum=0;
    for i=1:size(Cx,1)
        for j=1:size(Cx,2)
            if Cx(i,j)~=0 && prob_x(k,i)~=0 && prob_x(k,j)
                sum=sum+(-Cx(i,j)*log(Cx(i,j)/(prob_x(k,i)*prob_x(k,j))));
            end
        end
    end
    Ix(k)=sum;
end
