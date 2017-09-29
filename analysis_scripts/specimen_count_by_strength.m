function specimen_count_by_strength(database)

fid = 1;
database_folder = '..';
load(fullfile(database_folder,sprintf('%s.mat',database)));

fc = [data(:).fc];
Fy = [data(:).Fy];

fprintf(fid,'\n');
fprintf(fid,'                | Normal Strength Concrete |  High Strength Concrete  |  Total |\n');
fprintf(fid,'                |       f''c <= 10 ksi      |     f''c > 10 ksi         |        |\n');
fprintf(fid,'----------------+--------------------------+--------------------------+--------|\n');
fprintf(fid,'Normal Strength |                          |                          |        |\n');
fprintf(fid,'     Steel      |          %4i            |          %4i            |  %4i  |\n',sum(fc<=10 & Fy <= 75),sum(fc>10 & Fy <= 75),sum(Fy <= 75));
fprintf(fid,'  Fy <= 75 ksi  |                          |                          |        |\n');
fprintf(fid,'----------------+--------------------------+--------------------------+--------|\n');
fprintf(fid,' High Strength  |                          |                          |        |\n');
fprintf(fid,'     Steel      |          %4i            |          %4i            |  %4i  |\n',sum(fc<=10 & Fy > 75),sum(fc>10 & Fy > 75),sum(Fy > 75));
fprintf(fid,'  Fy > 75 ksi   |                          |                          |        |\n');
fprintf(fid,'----------------+--------------------------+--------------------------+--------|\n');
fprintf(fid,'     Total      |          %4i            |          %4i            |  %4i  |\n',sum(fc<=10),sum(fc>10),length(fc));
fprintf(fid,'\n');

end
