% input var a denotes which HF channel to take.

function [HF, HF1] = full_signal(loc,param,a,bad)
k = param.velmex.XNStep*param.velmex.YNStep;
if param.velmex.FastAxis == 'X'
    fL = param.velmex.XNStep; % gets fast direction scan points
    sL = param.velmex.YNStep; %gets slow direction scan points
else
    fL = param.velmex.YNStep;
    sL = param.velmex.XNStep;
end


%b = waitbar(0);
if bad
    sL = sL-1;
end
%     for i = 1:fL
%         for j = 1:sL
%             if read
%             [~,HF1{i,j}] = read_ucsdi_data(loc,(i-1)*sL+j); % Gets data sequentially
%         end
%         %waitbar(i/(fL),b,'Creating 4D array');
%         % fprintf('.');
%         multiWaitbar('Creating 4D Array',i/fL);
%     end

for j = 1:sL
    for i = 1:fL
        [~,HF1{i,j}] = read_ucsdi_data(loc,(i)+(fL*(j-1))); % Gets data sequentially
       % disp((i)+(fL*(j-1)));
    end
    %waitbar(i/(fL),b,'Creating 4D array');
    % fprintf('.');
    multiWaitbar('Creating 4D Array',i/sL);
end

% for i = 1:size(HF1,1)
%     for j = 1:size(HF1,2)
%     y(i,:,:) = HF1{i,1};
%     y2(j,:,:) = HF1{75,j};
%     end
% end

  

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



