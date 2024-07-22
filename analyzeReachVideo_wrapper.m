function analyzeReachVideo_wrapper(doFunction,vars)

switch doFunction 
    case 'extractEventsFromMovie'
        videoFile=vars.videoFile;
        endofVfname=vars.endofVfname;
        endofDir=vars.endofDir;
        zoneVals=vars.zoneVals;
        if vars.subtractExternalCue==true
            %pause;
            result_signal=removeExternalCue(zoneVals.cueZone,zoneVals.pelletZone,zoneVals.LEDZone,true);
            zoneVals.pelletZone(result_signal==0)=0;
            zoneVals.pelletZone(result_signal==1)=1;
        end
        [out,zoneVals,reaches,pellets,eat,paw,fidget,settings]=extractEventsFromMovie([videoFile(1:endofVfname(end)-1) '_zones.mat'],videoFile,zoneVals);
        
        savehandles=reformatAutoClassifyOutput(out,zoneVals,reaches,pellets,eat,paw,fidget,settings);
        save([videoFile(1:endofVfname(end)-1) '_savehandles.mat'],'savehandles');
    case 'getAlignment'
        savehandles=vars.savehandles;
        out=vars.out;
        videoFile=vars.videoFile;
        endofVfname=vars.endofVfname;
        fractionRange=vars.fractionRange;
        isInSecondHalf=vars.isInSecondHalf;
        fractionThroughArduino=vars.fractionThroughArduino;
        tryinc=vars.tryinc;
        try_scale1=vars.try_scale1;
        try_scale2=vars.try_scale2;
        
        % Discard end of video
        savehandles=discardLastNFrames(savehandles);
        
        % Align Arduino output data and data from video file
        if isfield(savehandles,'discardLastN')
            discardLastN=savehandles.discardLastN;
        else
            discardLastN=0;
        end
        settings=alignmentSettings(discardLastN,fractionRange,isInSecondHalf,fractionThroughArduino,[tryinc try_scale1 try_scale2]);
        aligned=getAlignment(out,30,savehandles,settings);
        save([videoFile(1:endofVfname(end)-1) '_aligned.mat'],'aligned');
        save([videoFile(1:endofVfname(end)-1) '_alignmentSettings.mat'],'settings');
    case 'getCue'
        minProm=vars.minProm;
        aligned2=vars.aligned2;
        videoFile=vars.videoFile;
        endofVfname=vars.endofVfname;
        [aligned,cleanup]=cleanUpCue_basedOnArduino(aligned2,minProm,minProm/2);
        save([videoFile(1:endofVfname(end)-1) '_cleanup_settings.mat'], 'cleanup');
        save([videoFile(1:endofVfname(end)-1) '_aligned.mat'],'aligned');
    case 'organizeData'
        warning off;
        videoFile=vars.videoFile;
        endofVfname=vars.endofVfname;
        savehandles=vars.savehandles;
        out=vars.out;
        aligned=vars.aligned;
        
        [status]=mkdir([videoFile(1:endofVfname(end)-1) '_processed_data']);
        finaldata=integrateSDoutWithReaches(savehandles,out,30,aligned,[videoFile(1:endofVfname(end)-1) '_processed_data']);
        
        % Check for chewing of pellet (this should take a while -- pellet is large)
        finaldata=checkForChewedPellet(finaldata);
        alignment=finaldata;
        save([videoFile(1:endofVfname(end)-1) '_processed_data' '/final_aligned_data.mat'],'alignment');
        
        % Plot results
        finaldata=alignment;
        if ~isfield(finaldata,'isGrooming')
            finaldata.isGrooming=zeros(size(finaldata.cueZone_onVoff));
        end
        tbt=plotCueTriggeredBehavior(finaldata,'cueZone_onVoff',0);
        save([videoFile(1:endofVfname(end)-1) '_processed_data/tbt.mat'],'tbt');
        settings=plotCueTriggered_settings();
        save([videoFile(1:endofVfname(end)-1) '_plottingSettings.mat'],'settings');
    otherwise
        error('Unrecognized value of doFunction argument to analyzeReachVideo_wrapper');
end
        
end

