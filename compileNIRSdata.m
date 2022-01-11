%Compiles the data for two subjects within a dyad and two conversations each. 
%Compiles your choice of channel rejection from preprocessing (i.e. ch_reject) 

function [deoxy3D,oxy3D]= compileNIRSdata(preprocess_dir,dataprefix,ch_reject,numScans)

%find all of the preprocessed folders
currdir=dir(strcat(preprocess_dir,filesep,dataprefix,'*'));

% loop to compile the data
    for i=1:length(currdir)
        dyad=currdir(i).name; %define dyad
        dyaddir=dir(strcat(preprocess_dir,filesep,dyad,filesep,dataprefix,'*'));

        for sj=1:length(dyaddir)
            subject=dyaddir(sj).name;
            subfolder=dir(strcat(currdir(1).folder,filesep,dyad,filesep,subject,filesep,dataprefix,'*')); 

            if ~isempty(subfolder) && sj==1
                subfiles=dir(strcat(subfolder(1).folder,filesep,subfolder(1).name,filesep,'*.mat'));
                load(strcat(subfiles(1).folder,filesep,subfiles(ch_reject).name)) 
                length_convo=length(z_oxy);
                z_deoxy1_1(1:length_convo,:,i)=z_deoxy(1:length_convo,:);
                z_oxy1_1(1:length_convo,:,i)=z_oxy(1:length_convo,:);

                subfiles=dir(strcat(subfolder(2).folder,filesep,subfolder(2).name,filesep,'*.mat'));
                load(strcat(subfiles(1).folder,filesep,subfiles(ch_reject).name))
                length_convo=length(z_oxy);
                z_deoxy1_2(1:length_convo,:,i)=z_deoxy(1:length_convo,:);
                z_oxy1_2(1:length_convo,:,i)=z_oxy(1:length_convo,:);
                
                if numScans>2
                    scansleft=length(subfolder)-2;
                    if scansleft>0
                        subfiles=dir(strcat(subfolder(3).folder,filesep,subfolder(3).name,filesep,'*.mat'));
                        load(strcat(subfiles(1).folder,filesep,subfiles(ch_reject).name)) 
                        length_convo=length(z_oxy);
                        z_deoxy1_3(1:length_convo,:,i)=z_deoxy(1:length_convo,:);
                        z_oxy1_3(1:length_convo,:,i)=z_oxy(1:length_convo,:);
                        
                        if scansleft>1
                            subfiles=dir(strcat(subfolder(4).folder,filesep,subfolder(4).name,filesep,'*.mat'));
                            load(strcat(subfiles(1).folder,filesep,subfiles(ch_reject).name)) 
                            length_convo=length(z_oxy);
                            z_deoxy1_4(1:length_convo,:,i)=z_deoxy(1:length_convo,:);
                            z_oxy1_4(1:length_convo,:,i)=z_oxy(1:length_convo,:);
                            
                            if scansleft>2
                                subfiles=dir(strcat(subfolder(5).folder,filesep,subfolder(5).name,filesep,'*.mat'));
                                load(strcat(subfiles(1).folder,filesep,subfiles(ch_reject).name)) 
                                length_convo=length(z_oxy);
                                z_deoxy1_5(1:length_convo,:,i)=z_deoxy(1:length_convo,:);
                                z_oxy1_5(1:length_convo,:,i)=z_oxy(1:length_convo,:);
                            end
                        end
                    end
                end

            elseif ~isempty(subfolder) && sj==2
                subfiles=dir(strcat(subfolder(1).folder,filesep,subfolder(1).name,filesep,'*.mat'));
                load(strcat(subfiles(1).folder,filesep,subfiles(ch_reject).name))
                length_convo=length(z_oxy);
                z_deoxy2_1(1:length_convo,:,i)=z_deoxy(1:length_convo,:);
                z_oxy2_1(1:length_convo,:,i)=z_oxy(1:length_convo,:);

                subfiles=dir(strcat(subfolder(2).folder,filesep,subfolder(2).name,filesep,'*.mat'));
                load(strcat(subfiles(1).folder,filesep,subfiles(ch_reject).name))
                length_convo=length(z_oxy);
                z_deoxy2_2(1:length_convo,:,i)=z_deoxy(1:length_convo,:);
                z_oxy2_2(1:length_convo,:,i)=z_oxy(1:length_convo,:);
                
                if numScans>2
                    scansleft=length(subfolder)-2;
                    if scansleft>0
                        subfiles=dir(strcat(subfolder(3).folder,filesep,subfolder(3).name,filesep,'*.mat'));
                        load(strcat(subfiles(1).folder,filesep,subfiles(ch_reject).name)) 
                        length_convo=length(z_oxy);
                        z_deoxy2_3(1:length_convo,:,i)=z_deoxy(1:length_convo,:);
                        z_oxy2_3(1:length_convo,:,i)=z_oxy(1:length_convo,:);

                        if scansleft>1
                            subfiles=dir(strcat(subfolder(4).folder,filesep,subfolder(4).name,filesep,'*.mat'));
                            load(strcat(subfiles(1).folder,filesep,subfiles(ch_reject).name)) 
                            length_convo=length(z_oxy);
                            z_deoxy2_4(1:length_convo,:,i)=z_deoxy(1:length_convo,:);
                            z_oxy2_4(1:length_convo,:,i)=z_oxy(1:length_convo,:);

                            if scansleft>2
                                subfiles=dir(strcat(subfolder(5).folder,filesep,subfolder(5).name,filesep,'*.mat'));
                                load(strcat(subfiles(1).folder,filesep,subfiles(ch_reject).name)) 
                                length_convo=length(z_oxy);
                                z_deoxy2_5(1:length_convo,:,i)=z_deoxy(1:length_convo,:);
                                z_oxy2_5(1:length_convo,:,i)=z_oxy(1:length_convo,:);
                            end
                        end
                    end
                end
            end
        end
    end
    deoxy3D(1).name='Scan 1';
    deoxy3D(1).sub1=z_deoxy1_1;
    deoxy3D(1).sub2=z_deoxy1_2;
    
    deoxy3D(2).name='Scan 2';
    deoxy3D(2).sub1=z_deoxy2_1;
    deoxy3D(2).sub2=z_deoxy2_2;
    
    oxy3D(1).name='Scan 1';
    oxy3D(1).sub1=z_oxy1_1;
    oxy3D(1).sub2=z_oxy1_2;
    
    oxy3D(2).name='Scan 2';
    oxy3D(2).sub1=z_deoxy2_1;
    oxy3D(2).sub2=z_oxy2_2;
    
    if numScans>2
        deoxy3D(3).name='Scan 3';
        deoxy3D(3).sub1=z_deoxy1_3;
        deoxy3D(3).sub2=z_deoxy2_3;

        oxy3D(3).name='Scan 3';
        oxy3D(3).sub1=z_oxy1_3;
        oxy3D(3).sub2=z_oxy2_3;
        
        if numScans>3
            deoxy3D(4).name='Scan 4';
            deoxy3D(4).sub1=z_deoxy1_4;
            deoxy3D(4).sub2=z_deoxy2_4;

            oxy3D(4).name='Scan 4';
            oxy3D(4).sub1=z_oxy1_4;
            oxy3D(4).sub2=z_oxy2_4;
            
            if numScans>4
                deoxy3D(5).name='Scan 5';
                deoxy3D(5).sub1=z_deoxy1_5;
                deoxy3D(5).sub2=z_deoxy2_5;

                oxy3D(5).name='Scan 5';
                oxy3D(5).sub1=z_oxy1_5;
                oxy3D(5).sub2=z_oxy2_5;
            end
        end
    end
end

