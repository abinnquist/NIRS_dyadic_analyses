%Compiles the data for two subjects within a dyad and two conversations
%each. Compiles the "nouncertainchannels" from preprocessing. 

function [z_deoxy1_1,z_oxy1_1,z_totaloxy1_1,z_deoxy1_2,z_oxy1_2,...
    z_totaloxy1_2,z_deoxy2_1,z_oxy2_1,z_totaloxy2_1,z_deoxy2_2,...
    z_oxy2_2,z_totaloxy2_2]= compileNIRSdata(preprocess_dir,dataprefix,ch_reject)

%find all of the preprocessed folders
currdir=dir(strcat(preprocess_dir,filesep,dataprefix,'*'));

% write a loop to compile the data
    for i=1:length(currdir)
        dyad=currdir(i).name; %define dyad
        dyaddir=dir(strcat(preprocess_dir,filesep,dyad,filesep,dataprefix,'*'));

        for sj=1:length(dyaddir)
            subject=dyaddir(sj).name;
            subfolder=dir(strcat(currdir(1).folder,filesep,dyad,filesep,subject,filesep,dataprefix,'*')); 

            if ~isempty(subfolder) && sj==1
                subfiles=dir(strcat(subfolder(1).folder,filesep,subfolder(1).name,filesep,'*.mat'));
                load(strcat(subfiles(1).folder,filesep,subfiles(ch_reject).name)) %Compile the "nouncertain channels"
                length_convo=length(z_oxy);
                z_deoxy1_1(1:length_convo,:,i)=z_deoxy(1:length_convo,:);
                z_oxy1_1(1:length_convo,:,i)=z_oxy(1:length_convo,:);
                z_totaloxy1_1(1:length_convo,:,i)=z_totaloxy(1:length_convo,:);

                subfiles=dir(strcat(subfolder(2).folder,filesep,subfolder(2).name,filesep,'*.mat'));
                load(strcat(subfiles(1).folder,filesep,subfiles(ch_reject).name))
                length_convo=length(z_oxy);
                z_deoxy1_2(1:length_convo,:,i)=z_deoxy(1:length_convo,:);
                z_oxy1_2(1:length_convo,:,i)=z_oxy(1:length_convo,:);
                z_totaloxy1_2(1:length_convo,:,i)=z_totaloxy(1:length_convo,:);

            elseif ~isempty(subfolder) && sj==2
                subfiles=dir(strcat(subfolder(1).folder,filesep,subfolder(1).name,filesep,'*.mat'));
                load(strcat(subfiles(1).folder,filesep,subfiles(ch_reject).name))
                length_convo=length(z_oxy);
                z_deoxy2_1(1:length_convo,:,i)=z_deoxy(1:length_convo,:);
                z_oxy2_1(1:length_convo,:,i)=z_oxy(1:length_convo,:);
                z_totaloxy2_1(1:length_convo,:,i)=z_totaloxy(1:length_convo,:);

                subfiles=dir(strcat(subfolder(2).folder,filesep,subfolder(2).name,filesep,'*.mat'));
                load(strcat(subfiles(1).folder,filesep,subfiles(ch_reject).name))
                length_convo=length(z_oxy);
                z_deoxy2_2(1:length_convo,:,i)=z_deoxy(1:length_convo,:);
                z_oxy2_2(1:length_convo,:,i)=z_oxy(1:length_convo,:);
                z_totaloxy2_2(1:length_convo,:,i)=z_totaloxy(1:length_convo,:);
            end
        end
    end
end
