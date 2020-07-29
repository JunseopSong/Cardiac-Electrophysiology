% Calculate unipolar electric potential from 2D simulation data
% Jun-Seop Song

clc; clear all;

%%
pivotX = 50;
pivotY = 50;
h = 0.1;  % cm
dx = 0.025;  % cm

sigma_i = 2.5;
sigma_e = 20;
constant = sigma_i/(2*pi*sigma_e);

startTime = 0;
DATA_Length = 200;
size = 200;
printTimeInterval = 10;

Vm = zeros(size, size);
V = zeros(DATA_Length, 1);
dV = zeros(size, size);

n1 = 128;
ia1 = max(2, pivotX-n1);
ib1 = min(size-1, pivotX+n1);
ja1 = max(2, pivotY-n1);
jb1 = min(size-1, pivotY+n1);

%%
for k = 1:DATA_Length
    fid = fopen(['vm' num2str(k*printTimeInterval + startTime) '.txt']);
    data_tmp = fread(fid, size*size, 'double');
    fclose(fid);
    Vm(:,:) = reshape(data_tmp, size, size);
    
    for i = ia1:ib1
        for j = ja1:jb1
            rho = sqrt((i-pivotX)^2 + (j-pivotY)^2) * dx;
            f_rho = (rho^2-2*h^2)/(rho^2+h^2)^(5/2);
            dV(i,j) = Vm(i,j)*f_rho;
        end
    end
    U_real(k) = constant * sum(sum(dV)) * (dx)^2;
    U_approximation(k) = -Vm(pivotX,pivotY);
    
    disp(k);
end


%%
figure;
subplot(2,1,1);
plot(printTimeInterval:printTimeInterval:DATA_Length*printTimeInterval, U_real);
subplot(2,1,2);
plot(printTimeInterval:printTimeInterval:DATA_Length*printTimeInterval, U_approximation);
