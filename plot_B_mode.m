% x and z are sample coord pairs for start and last
% t is the time sample or sample range

function plot_B_mode(HF,x,z,t)
for j = t
    for i = x(1):x(2)
    data{j}(i+1-x(1),:) = HF{i,1}(z(1):z(2),j);
    end
end
for i = t
    figure;
    imagesc(data{i}')
end
    

