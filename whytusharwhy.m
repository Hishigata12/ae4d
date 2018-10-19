function L = whytusharwhy(X)
A = table2array(X);
test = length(str2num(A(1,:)));
q = length(A);
for i = 1:q%test
    P(i,:) = str2num(A(i,:));
end
L = reshape(P,[test,test,test]);
% 
% for i = 1:81
%     imagesc(L(:,:,i))
%     drawnow
%     pause(0.1)
% end
