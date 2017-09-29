function specimen_count()

fid = 1;
database_folder = '..';
databases = {
    'CCFT_C+PBC'
    'RCFT_C+PBC'
     'SRC_C+PBC'
    'CCFT_Beams'
    'RCFT_Beams'
    'CCFT_Other'
    'RCFT_Other'};

% Count and print specimen totals
total = 0;
for i = 1:length(databases)
    load(fullfile(database_folder,sprintf('%s.mat',databases{i})));
    fprintf(fid,'%16s  %4i\n',databases{i},length(data));
    total = total + length(data);
end
fprintf(fid,'           Total  %4i\n',total);

end