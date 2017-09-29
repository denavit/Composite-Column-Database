function specimen_summary

[UniqueRefs1,UniqueRefCount1,UniqueRefCountStr1,UniqueAlpha1] = getData('CCFT','C+PBC');
[UniqueRefs2,UniqueRefCount2,UniqueRefCountStr2,UniqueAlpha2] = getData('RCFT','C+PBC');
[UniqueRefs3,UniqueRefCount3,UniqueRefCountStr3,UniqueAlpha3] = getData( 'SRC','C+PBC');
[UniqueRefs4,UniqueRefCount4,UniqueRefCountStr4,UniqueAlpha4] = getData('CCFT','Beams');
[UniqueRefs5,UniqueRefCount5,UniqueRefCountStr5,UniqueAlpha5] = getData('CCFT','Other');
[UniqueRefs6,UniqueRefCount6,UniqueRefCountStr6,UniqueAlpha6] = getData('RCFT','Beams');
[UniqueRefs7,UniqueRefCount7,UniqueRefCountStr7,UniqueAlpha7] = getData('RCFT','Other');

UniqueRefs          = horzcat(UniqueRefs1,UniqueRefs2,UniqueRefs3,...
    UniqueRefs4,UniqueRefs5,UniqueRefs6,UniqueRefs7);
UniqueRefCount      = horzcat(UniqueRefCount1,UniqueRefCount2,UniqueRefCount3,...
    UniqueRefCount4,UniqueRefCount5,UniqueRefCount6,UniqueRefCount7);
UniqueRefCountStr   = horzcat(UniqueRefCountStr1,UniqueRefCountStr2,UniqueRefCountStr3,...
    UniqueRefCountStr4,UniqueRefCountStr5,UniqueRefCountStr6,UniqueRefCountStr7);
UniqueAlpha         = horzcat(UniqueAlpha1,UniqueAlpha2,UniqueAlpha3,...
    UniqueAlpha4,UniqueAlpha5,UniqueAlpha6,UniqueAlpha7);

[UUniqueRefs,ia,ic] = unique(UniqueRefs);
UUniqueAlpha = UniqueAlpha(ia);
UUniqueRefCount = nan(1,length(UUniqueRefs));
UUniqueRefCountStr = cell(1,length(UUniqueRefs));
for i = 1:length(UUniqueRefs)
    ind = find(ic==i);
    if numel(ind) == 1
        UUniqueRefCount(i) = UniqueRefCount(ind);
        UUniqueRefCountStr{i} = UniqueRefCountStr{ind};
    else
        UUniqueRefCount(i) = sum(UniqueRefCount(ind));
        temp = sprintf('%i Total:',UUniqueRefCount(i));
        for j = 1:length(ind)
            temp = sprintf('%s\n\t%s',temp,UniqueRefCountStr{ind(j)});
        end
        UUniqueRefCountStr{i} = temp;
    end
end


% Rearrange by year then name
[~,ix] = sort(UUniqueAlpha);
UUniqueRefs          = UUniqueRefs(ix);
UUniqueRefCount      = UUniqueRefCount(ix);
UUniqueRefCountStr   = UUniqueRefCountStr(ix);



% Print
for i = 1:length(UUniqueRefs)  
    fprintf('%s\t%s\n',UUniqueRefs{i},UUniqueRefCountStr{i});
end

end




function [UniqueRefs,UniqueRefCount,UniqueRefCountStr,UniqueAlpha] = getData(sectionType,memberType)

database_folder = '..';
database_name = sprintf('%s_%s.mat',sectionType,memberType);
load(fullfile(database_folder,database_name))

numSpecimens = length(data);
Reference = cell(1,numSpecimens);
Alpha = cell(1,numSpecimens);

for i = 1:length(data)
    Reference{i} = sprintf('%s %s',data(i).Author,data(i).Year);
    Alpha{i}     = sprintf('%s %s',data(i).Year,data(i).Author);
end

% Find unique references and count
[UniqueRefs,ia,ic] = unique(Reference);
UniqueAlpha = Alpha(ia);
UniqueRefCount = nan(1,length(UniqueRefs));
UniqueRefCountStr = cell(1,length(UniqueRefs));
for i = 1:length(UniqueRefs)
    UniqueRefCount(i) = sum(ic==i);
    UniqueRefCountStr{i} = sprintf('%i %s %s',UniqueRefCount(i),sectionType,memberType);
end

end