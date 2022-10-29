function alignment=integrateSDoutWithReaches(reaches,out,moviefps,alignment,savedir)

% Will need to fix some file alignments by lining up distractor LED onsets and offsets
if isempty(alignment)
    alignment=getAlignment(out,moviefps,reaches);
    save([savedir '/alignment.mat'], 'alignment');
end

% Remove incomplete reach detections
if length(reaches.eatTime)<length(reaches.reachStarts)
    reaches.eatTime=[reaches.eatTime ones(1,length(reaches.reachStarts)-length(reaches.eatTime))*reaches.eatTime(end)];
end
isnotnaninds=~isnan(mean([reaches.reachStarts' reaches.pelletTime' reaches.eatTime'],2)) & (reaches.pelletTime'-reaches.reachStarts'>=0);
f=fieldnames(reaches);
l=length(reaches.reachStarts);
for i=1:length(f)
    if length(reaches.(f{i}))==l
        temp=reaches.(f{i});
        reaches.(f{i})=temp(isnotnaninds);
    end
end

% Flip pellet present because actually saved 1 if pellet MISSING
reaches.pelletPresent=reaches.pelletMissing==0;

% Add other events 
alignment=integrateReachTypes(reaches, alignment);

% Save data
save([savedir '/final_aligned_data.mat'],'alignment');

end
