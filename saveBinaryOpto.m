function tbt=saveBinaryOpto(tbt,settings)

% Save a binary opto field
optoZoneInd=nan;
for i=1:length(settings.plotevents)
    if strcmp(settings.plotevents{i},'optoZone')
        optoZoneInd=i;
        break
    end
end

if isnan(optoZoneInd)
    return
end

tbt.optoZone_binary=tbt.optoZone>settings.eventThresh{optoZoneInd};