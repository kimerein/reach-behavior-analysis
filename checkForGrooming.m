function eat=checkForGrooming(eat,settings)

movieframeinds=settings.discardFirstNFrames:settings.discardFirstNFrames+length(eat.isChewing)-1;

eat.isChewing_backup=eat.isChewing;

chew=eat.isChewing;
[labeledVector,numRegions]=bwlabel(chew);
% Check each stretch for grooming
disp('For each stretch of movie, is mouse grooming? Enter "yes", "no" or "both".');
isGrooming=zeros(size(chew));
i=1;
while i<=numRegions
    if any(isnan(chew(labeledVector==i)))
        % nans
    else
        % is this stretch grooming?
        fi=find(labeledVector==i,1,'first');
        li=find(labeledVector==i,1,'last');
        s=input(['Movie frame inds ' num2str(movieframeinds(fi)) ' to ' num2str(movieframeinds(li)) '. Grooming? '],'s');
        switch s
            case 'yes'
                % is grooming
                isGrooming(fi:li)=1;
            case 'no'
                % is not grooming
                isGrooming(fi:li)=0;
            case 'both'
                % includes both
                % have user separate
                sepf=input('Enter frame separating grooming from eating: ');
                if ~ismember(sepf,movieframeinds)
                    error('Entered movie frame is not contained in this video segment');
                else
                    [~,mf]=nanmin(abs(movieframeinds-sepf));
                    chew(mf)=0;
                    [labeledVector,numRegions]=bwlabel(chew);
                    continue % without incrementing i
                end
            otherwise
                error('Unrecognized user input. Input should be "yes", "no" or "both".');
        end
    end
    i=i+1;
end
isGrooming(isnan(chew))=nan;
eat.isGrooming=isGrooming;

if settings.removeGroomingFromEating==1
    % if is a grooming stretch, is not an eating stretch
    eat.isChewing(isGrooming==1)=0;
end
