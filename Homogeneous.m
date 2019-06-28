function [RegR] = Homogeneous(R)

yb = size(R,1);

Rsort(1,:) = R(1,:);

for i = 2:yb
    a = R(i,:);
    a = [a(i:yb),a(1:i-1)];
    Rsort(i,:) = a;
end
Mean = mean(Rsort,1);

% Make a circulant matrix with the mean correlation function
[row,column] = size(Mean);
if column ~= 1
    Mean = Mean';
    row = column;
end

RsortC = Mean;
for i = 1:row-1
    Mean = Mean([end 1:end-1]);
    RsortC =[RsortC Mean];
end

RegR = RsortC;