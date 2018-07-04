% input var a denotes which HF channel to take.

function [HF, HF1] = full_signal(loc,param,a)
k = param.velmex.XNStep*param.velmex.YNStep;
if param.velmex.FastAxis == 'X'
    fL = param.velmex.XNStep; % gets fast direction scan points
    sL = param.velmex.YNStep; %gets slow direction scan points
else
    fL = param.velmex.YNStep;
    sL = param.velmex.XNStep;
end


%b = waitbar(0);
for i = 1:fL
    for j = 1:sL
        [~,HF1{i,j}] = read_ucsdi_data(loc,(i-1)*sL+j); % Gets data sequentially
    end
    %waitbar(i/(fL),b,'Creating 4D array');
   % fprintf('.');
    multiWaitbar('Creating 4D Array',i/fL);
end

  

if length(size(HF1{1})) > 2
    if exist('a','var')
        for i = 1:size(HF1,1)
            for j = 1:size(HF1,2)
                HF1{i,j} = HF1{i,j}(:,:,a);
            end
        end
    end
end

HF=zeros(size(HF1,1),size(HF1,2),size(HF1{1},1),size(HF1{1},2));
% for i = 1:fL
%     for j = 1:sL
%        % HF((i-1)*sL+j,:,:) = HF1{i,j}; %converts cells to pseudo 4-D array
%        HF(i,j,:,:) = HF1{i,j}; %Converts cell array to double
%     end
%     waitbar(.5+i/(fL*2),b,'Finalizing 4D array construction');
% end

%HF = reshape(HF,fL,sL,param.daq.HFdaq.pts,param.daq.HFdaq.NoBurstTriggers); %converts pseudo 4-D to 4-D

if param.velmex.FastAxis == 'Y'
    HF = permute(HF,[2 1 3 4]); %rearrange Y and X if data was taken that way
end
%delete(b);
end



