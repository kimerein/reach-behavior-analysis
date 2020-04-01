function checkForControl(expt_dir)

% Need to check whether each experiment in this directory was non-standard in some way
% e.g.,
% Omitted cue
% Opto only (no cue)

ls=dir(expt_dir);
for i=1:length(ls)
    thisname=ls(i).name;
    thisisdir=ls(i).isdir;
    if (~isempty(regexp(thisname,'2017')) || ~isempty(regexp(thisname,'2018')) || ~isempty(regexp(thisname,'2019'))) && thisisdir==1
        resultsFiles=dir([expt_dir '\' thisname '\O2 output']);
        for j=1:length(resultsFiles)
            thisfilename=resultsFiles(j).name;
            thisfileisdir=resultsFiles(j).isdir;
            if ~isempty(regexp(thisfilename,'parsedOutput'))
                parsedStarts=regexp(thisfilename,'parsedOutput');
                getCurrentMovieName=thisfilename(1:parsedStarts-2);
                a=load([expt_dir '\' thisname '\O2 output\' thisfilename]);
                out=a.out;
                [wasControl,details]=checkForNoCueControl(out);
                if wasControl==1
                    % this was a no cue control
                    save([expt_dir '\' thisname '\O2 output\' getCurrentMovieName '_IS_NO_CUE_CONTROL.mat'],'details');
                    if exist([expt_dir '\' thisname '\O2 output\' getCurrentMovieName '_processed_data'], 'dir')==7
                        save([expt_dir '\' thisname '\O2 output\' getCurrentMovieName '_processed_data\IS_NO_CUE_CONTROL.mat'],'details');
                    end
                else
                    % was not no cue control
                    save([expt_dir '\' thisname '\O2 output\' getCurrentMovieName '_NOT_NO_CUE_CONTROL.mat'],'details');
                    if exist([expt_dir '\' thisname '\O2 output\' getCurrentMovieName '_processed_data'], 'dir')==7
                        save([expt_dir '\' thisname '\O2 output\' getCurrentMovieName '_processed_data\REGULAR.mat'],'details');
                    end
                end
                [wasControl,details]=checkForRedOptoOnly(out);
                if wasControl==1
                    % this was a red opto control
                    save([expt_dir '\' thisname '\O2 output\' getCurrentMovieName '_IS_RED_OPTO_CONTROL.mat'],'details');
                    if exist([expt_dir '\' thisname '\O2 output\' getCurrentMovieName '_processed_data'], 'dir')==7
                        save([expt_dir '\' thisname '\O2 output\' getCurrentMovieName '_processed_data\IS_RED_OPTO_CONTROL.mat'],'details');
                    end
                else
                    % was not red opto control
                    save([expt_dir '\' thisname '\O2 output\' getCurrentMovieName '_NOT_RED_OPTO_CONTROL.mat'],'details');
                    if exist([expt_dir '\' thisname '\O2 output\' getCurrentMovieName '_processed_data'], 'dir')==7
                        save([expt_dir '\' thisname '\O2 output\' getCurrentMovieName '_processed_data\REGULAR.mat'],'details');
                    end
                end
            end
        end
    end
end