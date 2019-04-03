%%acquisition of the needed data from the motion capture.
clearvars
myFolder = 'D:\stage\SO_entropy_based_motion_extraction_algo\fichiers c3d\cleaned\temp'; %define the folder where the files to compute are
files=get_file_names_c3d(myFolder); %Get the names of the files inside the folder
g=1;
for n=20:20:60
    g=g+3;
    for k=1:size(files,1)
        clearvars -except files k marker_set SAVE_DATA save_file_name first_frame g %clear all the variables except the one needed toexecute the loop over all the files
        file_name = files(k).name;  %get the name of the file that will be treated in this iteration
        acq = btkReadAcquisition(file_name); %read the file
        labels = fieldnames(btkGetMarkers(acq));   %get the labels of the markers
        markers_values = btkGetMarkersValues(acq); %give an array filled with the coordinates of the marker (x,y,z of first marker, x,y,z of second, and so on)
        events=btkGetEvents(acq);   %get the events stored in the c3d file (rigth foot strike...)
        first_frame(k)=btkGetFirstFrame(acq);   %Get the number of the first frame to correct the delay between the events detected by the algo and the event directly in the files
        
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
        n = 500; %number of discretizations levels, this setting may have importance on the results
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
        discretization_y=linspace(min_y,max_y,n);
        discretization_z=linspace(min_z,max_z,n);
        for i=1:size(displacement_x,1)    %will generate an histogram for each frame for the x,y and z coordinates
            histogram_x(i,:)=histcounts(displacement_x(i,:),discretization_x);
            histogram_y(i,:)=histcounts(displacement_y(i,:),discretization_y);
            histogram_z(i,:)=histcounts(displacement_z(i,:),discretization_z);  
        end
        
        %computing the probabilities
        for i=1:size(histogram_x)
            for j=1:size(histogram_x,2)
                prob_x(i,j)=histogram_x(i,j)/(size(markers_values,2)/3);
                prob_y(i,j)=histogram_y(i,j)/(size(markers_values,2)/3);
                prob_z(i,j)=histogram_z(i,j)/(size(markers_values,2)/3);
            end
        end         
        
        %filling the C(r,s) matrix
        %with the same global range as the histograms
        %create the matrices for each marker
        [matrix_r_s_global_x,matrix_r_s_x,Cx,proba_r_s_x]=compute_C(histogram_x,displacement_x,discretization_x,labels);
        [matrix_r_s_global_y,matrix_r_s_y,Cy,proba_r_s_y]=compute_C(histogram_y,displacement_y,discretization_y,labels);
        [matrix_r_s_global_z,matrix_r_s_z,Cz,proba_r_s_z]=compute_C(histogram_z,displacement_z,discretization_z,labels);
        
        %compute the mutual information for each dimension and sum it to get the
        %total mutual information I
        Ix=mutual_info(displacement,Cx,prob_x);
        Iy=mutual_info(displacement,Cy,prob_y);
        Iz=mutual_info(displacement,Cz,prob_z);
        I=Ix+Iy+Iz;
        
        subplot(5,3,g);
        plot(I);
        title('I');
        subplot(5,3,g+1);
        plot(Ix);
        title('Ix');
        subplot(5,3,g+2);
        plot(Iz);
        title('Iz');
        
        
        file_name=erase(file_name,'.c3d');
        save_file_name(k)=string(file_name);
        SAVE_DATA.(file_name).Ix=Ix;
        SAVE_DATA.(file_name).Iy=Iy;
        SAVE_DATA.(file_name).Iz=Iz;
        SAVE_DATA.(file_name).I=I;
        if isa(events,'struct')==1
            events=events.General_impact;
        end
        SAVE_DATA.(file_name).forceplate=events;
        
    end
end
%{
%set=["all","legs","shank","foot"];
window_size=20;
threshold=0;
clear comparison
comparison(1,1)=string('window size');
comparison(1,2)=window_size;
comparison(1,3)=string('threshold');
comparison(1,4)=threshold;
b=2;
for j=1:size(save_file_name,2)
    clear keyposes_temp keyposes_frames I_localized I_trimmed
    name=char(save_file_name(1,j));
    I=SAVE_DATA.(name).I;
    threshold=0;
    [keyposes_temp,I_localized,I_trimmed]=keyposes_detection(I,window_size,threshold);
    while size(keyposes_temp,2)>3
        [keyposes_temp,I_localized,I_trimmed]=keyposes_detection(I,window_size,threshold);
        threshold=threshold-1;
    end
    
    keyposes_frames(1:size(keyposes_temp,2))=keyposes_temp;
    SAVE_RESULTS.(name).keyposes_frames=keyposes_frames;
    SAVE_RESULTS.(name).I_localized=I_localized;
    SAVE_RESULTS.(name).I_trimmed=I_trimmed;
    %get the timing for the i_trimmed
    frames_global=keyposes_frames+first_frame(j);
    time_cine=frames_global/200;
    SAVE_RESULTS.(name).time_cine=time_cine;
    comparison(b,1)=string(name);
    b=b+1;
    comparison(b,1)='events: forceplate';
    b=b+1;
    comparison(b,1:size(SAVE_DATA.(name).forceplate,2))=SAVE_DATA.(name).forceplate;
    b=b+1;
    comparison(b,1)=('events: algo');
    b=b+1;
    comparison(b,1:size(time_cine,2))=time_cine;
    b=b+2;
    time_cine=[];
end
comparison=fillmissing(comparison,'constant',string(' '));


%plot(I_trimmed);
%scatter_3D_plot(markers_values,keyposes_frames)
%legend()

%}
%[frame_contact,frame_off,frame_contact2,frame_off2]=detect_event_forceplate(acq);
