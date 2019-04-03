%%acquisition of the needed data from the motion capture.
clearvars
myFolder = 'D:\stage\Sofamehack2019\all_files'; %define the folder where the files to compute are
files=get_file_names_c3d(myFolder); %Get the names of the files inside the folder
for l=1:4   %execute the whole treatment four times with different marker sets
    if l==1
        marker_set='all';   %takes in consideration all the markers
    else if l==2
            marker_set='legs';  %takes in consideration the markers of the leg
        else if l==3
                marker_set='shank'; %takes in consideration the markers of the shank
            else
                marker_set='foot';  %takes in consideration the markers of the foot
            end
        end
    end
    for k=1:size(files,1)
        clearvars -except files k marker_set SAVE_DATA save_file_name   %clear all the variables except the one needed toexecute the loop over all the files
        file_name = files(k).name;  %get the name of the file that will be treated in this iteration
        acq = btkReadAcquisition(file_name); %read the file        
        labels = fieldnames(btkGetMarkers(acq));   %get the labels of the markers
        markers_values = btkGetMarkersValues(acq); %give an array filled with the coordinates of the marker (x,y,z of first marker, x,y,z of second, and so on)
        events=btkGetEvents(acq);   %get the events stored in the c3d file (rigth foot strike...)
        if strcmp(marker_set,'legs')==1 %the data is organized in such a way that deleting a part of the rigth side of the array is enough to keep the markers of interest
            markers_values(:,1:17*3)=[];    %both side are considered, the markers are: ASI PSI THI KNE TIB ANK HEE TOE
        else if strcmp(marker_set,'shank')==1
                markers_values(:,1:23*3)=[];     %both side are considered, the markers are: KNE TIB ANK HEE TOE
            else if strcmp(marker_set,'foot')==1
                    markers_values(:,1:29*3)=[];    %both side are considered, the markers are: HEE TOE
                    for i=1:3
                        for j=1:size(markers_values,1)
                            temp_markers_values(j,i)=(markers_values(j,i)+markers_values(j,i+9))/2; %the midpoint of the markers of the foot is considered (right)
                        end
                    end
                    for i=4:6
                        for j=1:size(markers_values,1)
                            temp_markers_values(j,i)=(markers_values(j,i)+markers_values(j,i+3))/2; %the midpoint of the markers of the foot is considered (left)
                        end
                    end
                    markers_values=[];  %clear the whole array to copy the temporary array into it 
                    markers_values=temp_markers_values; %copy remaining data into the cleaned array
                end
            end
        end
        
        %%calculation of the displacement vector
        displacement=zeros(size(markers_values,1)-1,size(markers_values,2));    %preallocation
        for i=1:(size(markers_values,1)-1) %-1 otherwise it will crash trying to reach the last frame + 1            
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
        n = 257; %number of discretizations levels, this setting may have importance on the results
        maximum=max(displacement,[],'all');  %gives the max values for all the column hence for all the x,y,z coordinates respectively
        minimum=min(displacement,[],'all');  %same for min
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
        discretization_y=linspace(min_y,max_y,n);
        discretization_z=linspace(min_z,max_z,n);
        discretization_global=linspace(minimum,maximum,n);
        for i=1:size(displacement_x,1)    %will generate an histogram for each frame for the x,y and z coordinates
            histogram_x(i,:)=histcounts(displacement_x(i,:),discretization_x);
            histogram_y(i,:)=histcounts(displacement_y(i,:),discretization_y);
            histogram_z(i,:)=histcounts(displacement_z(i,:),discretization_z);
            histogram_max_min(i,:)=histcounts(displacement(i,:),discretization_global);
        end
        
        
        %computing the probabilities
        for i=1:size(histogram_x)
            for j=1:size(histogram_x,2)
                prob_x(i,j)=histogram_x(i,j)/(size(markers_values,2)/3);
                prob_y(i,j)=histogram_y(i,j)/(size(markers_values,2)/3);
                prob_z(i,j)=histogram_z(i,j)/(size(markers_values,2)/3);
                prob_global(i,j)=histogram_max_min(i,j)/(size(markers_values,2)/3);
            end
        end     
        
        
        %filling the C(r,s) matrix
        %with the same global range as the histograms
        %create the matrices for each marker
        [matrix_r_s_global_x,matrix_r_s_x,Cx,proba_r_s_x]=compute_C(histogram_x,displacement_x,discretization_x,labels);
        [matrix_r_s_global_y,matrix_r_s_y,Cy,proba_r_s_y]=compute_C(histogram_y,displacement_y,discretization_y,labels);
        [matrix_r_s_global_z,matrix_r_s_z,Cz,proba_r_s_z]=compute_C(histogram_z,displacement_z,discretization_z,labels);
        [matrix_r_s_global_max_min,matrix_r_s_max_min,C_max_min,proba_r_s_max_min]=compute_C(histogram_max_min,displacement,discretization_global,labels);
        %compute the mutual information for each dimension and sum it to get the
        %total mutual information I
        Ix=mutual_info(displacement,Cx,prob_x);
        Iy=mutual_info(displacement,Cy,prob_y);
        Iz=mutual_info(displacement,Cz,prob_z);
        I=Ix+Iy+Iz;
        I_max_min=mutual_info(displacement,C_max_min,prob_global);
        plot(I,'marker','o');
        hold on
        plot(Ix,'marker','+');
        hold on
        plot(Iy,'marker','*');
        hold on
        plot(Iz,'marker','x');
        hold on
        plot(I_max_min,'marker','s');
        legend('I','Ix','Iy','Iz','I_max_min');
        
        hold off
        file_name=erase(file_name,'.c3d');
        save_file_name(k)=string(file_name);
        SAVE_DATA.(marker_set).(file_name).Ix=Ix;
        SAVE_DATA.(marker_set).(file_name).Iy=Iy;
        SAVE_DATA.(marker_set).(file_name).Iz=Iz;
        SAVE_DATA.(marker_set).(file_name).I=I;
        
        %if events.Right_Foot_Strike_GS
        %{
        SAVE_DATA.(marker_set).(file_name).RFS=events.Right_Foot_Strike_GS;
        SAVE_DATA.(marker_set).(file_name).LFS=events.Left_Foot_Strike_GS;
        SAVE_DATA.(marker_set).(file_name).RFO=events.Right_Foot_Off_GS;
        SAVE_DATA.(marker_set).(file_name).LFO=events.Right_Foot_Off_GS;
        %}
    end
end
set=["all","legs","shank","foot"];
window_size=20;
threshold=0;
for i=1:4
    temp_set=set(i);
    for j=1:size(save_file_name,2)
        name=char(save_file_name(1,j));
        I=SAVE_DATA.(temp_set).(name).I;
        [keyposes_temp,I_localized,I_trimmed]=keyposes_detection(I,window_size,threshold);
        keyposes_frames(1:size(keyposes_temp,2))=keyposes_temp;
        SAVE_RESULTS.(temp_set).(name).keyposes_frames=keyposes_frames;
        SAVE_RESULTS.(temp_set).(name).I_localized=I_localized;
        SAVE_RESULTS.(temp_set).(name).I_trimmed=I_trimmed;
        %get the timing for the i_trimmed
        first_frame=btkGetFirstFrame(acq);
        frames_global=keyposes_frames+first_frame;
        time_cine=frames_global/200;
    end
end

%plot(I_trimmed);
%scatter_3D_plot(markers_values,keyposes_frames)
%legend()


%[frame_contact,frame_off,frame_contact2,frame_off2]=detect_event_forceplate(acq);
