function dbs_procData=checkForControl(varargin)

if length(varargin)==1
    expt_dir=varargin{1};
    dbs_procData=[];
elseif length(varargin)==2
    expt_dir=varargin{1};
    dbs_procData=varargin{2};
end

% Need to check whether each experiment in this directory was non-standard in some way
% e.g.,
% Omitted cue
% Opto only (no cue)

if ~iscell(expt_dir)
    ls=dir(expt_dir);
else
    ls=expt_dir;
end
areControls=zeros(1,length(ls));
backup_expt=expt_dir;
disp('Counting to ');
disp(length(ls));
for i=1:length(ls)
    disp(i);
    if ~iscell(backup_expt)
        thisname=ls(i).name;
        thisisdir=ls(i).isdir;
        totName=[expt_dir '\' thisname];
    else
        currdir=ls{i};
        temp=regexp(currdir,'\');
        thisname=currdir(temp(end)+1:end);
        thisisdir=isempty(regexp(thisname,'\.','ONCE'));
        totName=currdir;
        rar=regexp(currdir,'\');
        expt_dir=currdir(1:rar(end)-1);
    end
    oi=regexp(totName,'\');
    oneUpDir=totName(1:oi(end));
    if ~isempty(regexp(thisname,'processed_data','ONCE')) && thisisdir==1
    %if (~isempty(regexp(thisname,'2017')) || ~isempty(regexp(thisname,'2018')) || ~isempty(regexp(thisname,'2019'))) && thisisdir==1
        resultsFiles=dir(oneUpDir); 
        for j=1:length(resultsFiles)
            thisfilename=resultsFiles(j).name;
            thisfileisdir=resultsFiles(j).isdir;
            if ~isempty(regexp(thisfilename,'parsedOutput','ONCE'))
                parsedStarts=regexp(thisfilename,'parsedOutput');
                getCurrentMovieName=thisfilename(1:parsedStarts-2);
                a=load([expt_dir '\' thisfilename]);
                out=a.out;
                [wasControl,details]=checkForNoCueControl(out);
                if wasControl==1
                    % this was a no cue control
                    save([expt_dir  '\' getCurrentMovieName '_IS_NO_CUE_CONTROL.mat'],'details');
                    areControls(i)=1;
                    if exist([expt_dir  '\' getCurrentMovieName '_processed_data'], 'dir')==7
                        save([expt_dir  '\' getCurrentMovieName '_processed_data\IS_NO_CUE_CONTROL.mat'],'details');
                    end
                else
                    % was not no cue control
                    save([expt_dir  '\' getCurrentMovieName '_NOT_NO_CUE_CONTROL.mat'],'details');
                    if exist([expt_dir  '\' getCurrentMovieName '_processed_data'], 'dir')==7
                        save([expt_dir  '\' getCurrentMovieName '_processed_data\REGULAR.mat'],'details');
                    end
                end
                [wasControl,details]=checkForRedOptoOnly(out);
                if wasControl==1
                    % this was a red opto control
                    save([expt_dir  '\' getCurrentMovieName '_IS_RED_OPTO_CONTROL.mat'],'details');
                    areControls(i)=1;
                    if exist([expt_dir  '\' getCurrentMovieName '_processed_data'], 'dir')==7
                        save([expt_dir  '\' getCurrentMovieName '_processed_data\IS_RED_OPTO_CONTROL.mat'],'details');
                    end
                else
                    % was not red opto control
                    save([expt_dir  '\' getCurrentMovieName '_NOT_RED_OPTO_CONTROL.mat'],'details');
                    if exist([expt_dir  '\' getCurrentMovieName '_processed_data'], 'dir')==7
                        save([expt_dir  '\' getCurrentMovieName '_processed_data\REGULAR.mat'],'details');
                    end
                end
            end
        end
    end
end
if ~isempty(dbs_procData)
    dbs_procData.areControls=areControls;
end