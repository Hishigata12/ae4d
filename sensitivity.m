function sens = sensitivity(q,L,r)

%Creates range of HF data to look at
H  = q(r(1):r(end),:);

%Maps a point on the AE reconstruction to the LF data
c = linspace(1,length(L),size(H,2));
cmap = ceil(c);
L2 = L(cmap);

for i = 1:length(r)
sline{i} = abs(H(i,:))./abs(L2');
end

% Normalize values
Ln = L2/max(L2);
for i = 1:size(H,1)
Hn{i} = H(i,:)/max(H(i,:));
end

g = 6;


