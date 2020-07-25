clear all;

%%
utilized1 = load('utilized1.csv'); utilized2 = load('utilized2.csv'); utilized3 = load('utilized3.csv');
peaktopeak1 = load('peaktopeak1.csv'); peaktopeak2 = load('peaktopeak2.csv'); peaktopeak3 = load('peaktopeak3.csv');
xyz1 = load('xyz1.csv'); xyz2 = load('xyz2.csv'); xyz3 = load('xyz3.csv');

utilized = [utilized1 utilized2 utilized3];
peaktopeak = [peaktopeak1 peaktopeak2 peaktopeak3];
xyz = [xyz1 xyz2 xyz3];

sampleNum = length(utilized);
isExclude = zeros(sampleNum,1);

fid = fopen('Exclude.txt', 'w');
for i = 1:sampleNum
    if utilized(i) == 0
        isExclude(i) = 1;
    elseif peaktopeak(i) == 0 % Scar
        isExclude(i) = 1;
    end
    
    fprintf(fid, '%d\n', isExclude(i));
end
fclose(fid);