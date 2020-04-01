function summary=convertDatetimeToString(summary)

for i=1:size(summary.sess_datetime,1)
    summary.sess_datetime{i}=datestr(summary.sess_datetime{i});
    summary.trial_starts{i}=datestr(summary.trial_starts{i});
    summary.trial_ends{i}=datestr(summary.trial_ends{i});
end