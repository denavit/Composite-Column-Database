function test_to_predicted_ratio_by_reference(ttp_type,database)

fid = 1;
database_folder = '..';
load(fullfile(database_folder,sprintf('%s.mat',database)));

% Get references
refs = {data(:).Reference};
unique_refs = unique(refs);

% Get test to predicted ratio
switch ttp_type
    case 'AISC2016'
        ttp = [data(:).AISC2016_test_to_predicted];
    case 'PSD'
        ttp = [data(:).PSD_test_to_predicted];
    case 'ACDB'
        ttp = [data(:).ACDB_test_to_predicted];
    case 'Analysis_PfD'
        ttp = [data(:).Analysis_PfD_test_to_predicted];
    otherwise
        error('Unknown test to predicted type: %s',ttp_type)
end

% Compute and print results (formated to be copied to an excel spreadsheet)
fprintf(fid,'\t\tTest-to-predicted Ratio\t\t\t\n');
fprintf(fid,'Reference\tCount\tMin\tMax\tAvg\tStDev\n');
for i = 1:length(unique_refs)
    ind  = strcmp(refs,unique_refs{i});
    ittp = ttp(ind);
    fprintf(fid,'%s\t%i\t%g\t%g\t%g\t%g\n',unique_refs{i},sum(ind),...
        min(ittp),max(ittp),mean(ittp),std(ittp));
end

end