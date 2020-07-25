clear all;

%%
xyz1 = load('xyz1.csv'); ShEn1 = load('ShEn1.txt'); CFAE1 = load('CFAE1.csv')';
xyz2 = load('xyz2.csv'); ShEn2 = load('ShEn2.txt'); CFAE2 = load('CFAE1.csv')';
xyz3 = load('xyz3.csv'); ShEn3 = load('ShEn3.txt'); CFAE3 = load('CFAE1.csv')';

xyz = [xyz1 xyz2 xyz3];
ShEn_Raw = [ShEn1' ShEn2' ShEn3'];
CFAE_Raw = [CFAE1' CFAE2' CFAE3'];
clear xyz1 xyz2 xyz3 ShEn1 ShEn2 ShEn3 CFAE1 CFAE2 CFAE3;

sampleNum = length(xyz);

fid = fopen('ShEn_Raw_Data.txt', 'w');
for i = 1:sampleNum
    if CFAE_Raw(i) > 500
        CFAE_Raw(i) = NaN;
        ShEn_Raw(i) = NaN;
    end
    
    fprintf(fid, '%f %f %f %f\n', xyz(1,i), xyz(2,i), xyz(3,i), ShEn_Raw(i));
end
fclose(fid);


%%
node = load('CT_Node.txt');
element = load('CT_Element.txt');
nde = length(node); nel = length(element);
ShEn = zeros(nde, 1);
isProjectionNode = zeros(nde, 1);

%% Projection
for j = 1:sampleNum
    if isnan(ShEn_Raw(j))
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
    ShEn(minId) = ShEn_Raw(j);
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
    
    sum1 = 0; sum2 = 0;
    
    for j = 1:sampleNum
        if isnan(ShEn_Raw(j))
            continue;
        end
        
        d = sqrt((node(i,1)-xyz(1,j))^2 + (node(i,2)-xyz(2,j))^2 + (node(i,3)-xyz(3,j))^2);
        weight = (max(0, R-d)/(R*d))^2;
        sum1 = sum1 + weight*ShEn_Raw(j);
        sum2 = sum2 + weight;
    end
    
    ShEn(i) = sum1 / sum2;
end

%% Export
fid = fopen('ShEn_Map.plt', 'w');
fprintf(fid, 'VARIABLES = "X", "Y", "Z", "ShEn"\nZONE F=FEPOINT, ET=triangle, N=%d , E=%d\n', nde, nel);
for i = 1:nde
    if isnan(ShEn(i))
        fprintf(fid, '%f %f %f 0\n', node(i,1), node(i,2), node(i,3));
    else
        fprintf(fid, '%f %f %f %f\n', node(i,1), node(i,2), node(i,3), ShEn(i));
    end
end
for i = 1:nel
    fprintf(fid, '%d %d %d\n', element(i,1), element(i,2), element(i,3));
end
fclose(fid);