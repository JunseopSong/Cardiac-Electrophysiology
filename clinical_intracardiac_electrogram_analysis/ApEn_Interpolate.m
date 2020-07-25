clear all;

%%
xyz1 = load('xyz1.csv'); ApEn1 = load('ApEn1.txt');
xyz2 = load('xyz2.csv'); ApEn2 = load('ApEn2.txt');
xyz3 = load('xyz3.csv'); ApEn3 = load('ApEn3.txt');
xyz4 = load('xyz4.csv'); ApEn4 = load('ApEn4.txt');

xyz = [xyz1 xyz2 xyz3 xyz4];
ApEn_Raw = [ApEn1' ApEn2' ApEn3' ApEn4'];
clear xyz1 xyz2 xyz3 xyz4 ApEn1 ApEn2 ApEn3 ApEn4;

sampleNum = length(xyz);

isExclude = load('Exclude.txt');

fid = fopen('ApEn_Raw_Data.txt', 'w');
for i = 1:sampleNum
    if isExclude(i) == 1
        ApEn_Raw(i) = NaN;
    end
    
    fprintf(fid, '%f %f %f %f\n', xyz(1,i), xyz(2,i), xyz(3,i), ApEn_Raw(i));
end
fclose(fid);


%%
node = load('CT_Node.txt');
element = load('CT_Element.txt');
nde = length(node); nel = length(element);
ApEn = zeros(nde, 1);
isProjectionNode = zeros(nde, 1);

%% Projection
for j = 1:sampleNum
    if isnan(ApEn_Raw(j))
        continue;
    end
    
    minD = 1000000; minId = 1;
    
    for i = 1:nde
        d = (node(i,1)-xyz(1,j))^2 + (node(i,2)-xyz(2,j))^2 + (node(i,3)-xyz(3,j))^2;
        if d < minD
            minD = d;
            minId = i;
        end
    end
    
    isProjectionNode(minId) = 1;
    ApEn(minId) = ApEn_Raw(j);
    for i = 1:3
        xyz(i,j) = node(minId,i);
    end
end

%% Interpolate
weight = zeros(sampleNum, 1);
R = 10;

for i = 1:nde
    if isProjectionNode(i) == 1
        continue;
    end
    
    sum1 = 0; sum2 = 0; measureCnt = 0;
    
    for j = 1:sampleNum
        if isnan(ApEn_Raw(j))
            continue;
        end
        
        measureCnt = measureCnt + 1;
        
        d = sqrt((node(i,1)-xyz(1,j))^2 + (node(i,2)-xyz(2,j))^2 + (node(i,3)-xyz(3,j))^2);
        weight = (max(0, R-d)/(R*d))^2;
        sum1 = sum1 + weight*ApEn_Raw(j);
        sum2 = sum2 + weight;
    end
    
    if measureCnt <= 1
        ApEn(i) = NaN;
    else
        ApEn(i) = sum1 / sum2;
    end
end

%% Export
fid = fopen('ApEn_Map.plt', 'w');
fprintf(fid, 'VARIABLES = "X", "Y", "Z", "ApEn"\nZONE F=FEPOINT, ET=triangle, N=%d , E=%d\n', nde, nel);
for i = 1:nde
    if isnan(ApEn(i))
        fprintf(fid, '%f %f %f 0\n', node(i,1), node(i,2), node(i,3));
    else
        fprintf(fid, '%f %f %f %f\n', node(i,1), node(i,2), node(i,3), ApEn(i));
    end
end
for i = 1:nel
    fprintf(fid, '%d %d %d\n', element(i,1), element(i,2), element(i,3));
end
fclose(fid);