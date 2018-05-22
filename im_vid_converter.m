function movi = im_vid_converter(im,f,a)
[file path] = uigetfile(fullfile(pwd,'*.png*'));
[A,map,alpha] = imread([path file]); %AE file
for i = 1:f
    F(:,:,:,i) = A;
end

short = file(1:(end-4));
imwrite(A, [path file],'Alpha', ones(size(A,1),size(A,2)));
clear A map alpha
[A,map,alpha] = imread([path file]);

[file2 path2] = uigetfile(fullfile(pwd,'*png'));
B = imread([path2 file2]); %PE file
imwrite(A, [path2 file2],'Alpha', ones(size(A,1),size(A,2)));
clear B
B = imread([path2 file2]);
% Edit to splice in a second image
% for i = 2:2:f
%     F(:,:,:,i) = B;
% end

%Edit to add transparency and resave
short = file(1:(end-4)); short2 = file2(1:(end-4));
g = imshow(B);
hold 
imshow(A);
if exist('a','var') ==1
set(g,  'AlphaData', a);
saveas(g,[path short2 '_t.png']);
end
hold;
D = imread([path2 short2 '_t.png']);
hold;
C = imread([path short '_t.png']);
imshow(D)
hold
imshow(C)
%%
movi = immovie(F);
implay(F);
movie2avi(F,[path short '_m.avi'],'Compression','Cinepak')
V = VideoWriter(movi,'Uncompressed AVI');
open(V);
writeVideo(V,F);
close(V);
