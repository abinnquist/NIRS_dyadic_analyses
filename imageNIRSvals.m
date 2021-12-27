% Imaging correlated areas
function fileMade=imageNIRSvals(mniCoords)

        % Add SPM12 and xjview to use for imaging
    genPath=uigetdir('','Select the folder path for SPM and xjview');
    addpath(genpath(strcat(genPath,filesep,'xjview')))
    addpath(genpath(strcat(genPath,filesep,'spm12')))
    
    [dataValues, dataPath] = uigetfile('*.mat','Choose the values you want to image');
    data2image = strcat(dataPath,dataValues);
    data2image = struct2array(load(data2image));
    
    %Choose the dyad and conversation you wish to image
    convoQ = 'Do you want to image the first (1) disccusion or second (2) Conversation? \n';
    convo = input(convoQ);
    
    if convo==1
        convo2img(:,:)=data2image(:,:,1);
        conversation='Convo_1';
    else
        convo2img(:,:)=data2image(:,:,2);
        conversation='Convo_2';
    end
    
    dyadQ = 'What dyad would you like to image (1 to n)? \n';
    dyad2img = input(dyadQ);

    imgName=strcat('dyad_',num2str(dyad2img),'_',conversation,'.img');

    % Which conversation correlations you want to visualize
    convo_mask=convo2img(dyad2img,:)';
    
    % Make sure to change the name of the file
    nirs2img(imgName, mniCoords, convo_mask, 0,0,0)
    
    fileMade='HDR file for imaging saved to current directory';
end
