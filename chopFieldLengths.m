function summary=chopFieldLengths(summary)

% 32 is max field length
f=fieldnames(summary);
for i=1:length(f)
    currfield=f{i};
    if length(currfield)>31
        % need to chop
        summary.(currfield([1:14 end-14+1:end]))=summary.(currfield);
        summary=rmfield(summary,currfield);
    end
end