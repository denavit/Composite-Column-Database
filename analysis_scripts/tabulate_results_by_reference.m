clear all; close all; clc;

database = 'RCFT_C+PBC';
reference = 'Vrcelj & Uy 2002';
data_to_print = {'H','t','compactness','AISC2016_test_to_predicted'};
print_format = {'%.3f','%.3f','%s','%.3f'};

fid = 1;
database_folder = '..';
load(fullfile(database_folder,sprintf('%s.mat',database)));

ind = strcmp(reference,{data(:).Reference});


fprintf(fid,'Specimen');
for i = 1:length(data_to_print)
    fprintf(fid,'\t%s',data_to_print{i});
end
fprintf(fid,'\n');
for i = find(ind == 1)
    fprintf(fid,'%s',data(i).Specimen);
    for j = 1:length(data_to_print)
        fprintf(fid,'\t');
        fprintf(fid,print_format{j},data(i).(data_to_print{j}));
    end
    fprintf(fid,'\n');
end