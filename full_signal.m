% input var a denotes which HF channel to take.
% location is the filepath of the ae files
% param is the parameter file
% a designates the channel to process
% one: if 1 designates single element transducer
% new: if 1 designates that scan was taken using new AEScan

function [HF, HF1] = full_signal(loc,param,a)
k = param.velmex.XNStep*param.velmex.YNStep;
% if param.velmex.FastAxis == 'X'
    fL = param.velmex.XNStep; % gets fast direction scan points
    sL = param.velmex.YNStep; % gets slow direction scan points
% else
%     fL = param.velmex.YNStep;
%     sL = param.velmex.XNStep;
% end
one = param.post.onemhz;
new = param.post.new;
ind = param.post.ind;

%b = waitbar(0);
%     for i = 1:fL
%         for j = 1:sL
%             if read
%             [~,HF1{i,j}] = read_ucsdi_data(loc,(i-1)*sL+j); % Gets data sequentially
%         end
%         %waitbar(i/(fL),b,'Creating 4D array');
%         % fprintf('.');
%         multiWaitbar('Creating 4D Array',i/fL);
%     end
if one
    for j = 1:sL
        for i = 1:fL
            if new
                [~,HF1(i,j,:,:)] = Read_Data2(loc,(i-1)*sL+j,param); % Gets data sequentially
            else
                [~,HF1(i,j,:,:)] = read_ucsdi_data(loc,(i-1)*sL+j); % Gets data sequentially
            end
            % disp((i)+(fL*(j-1)));
            multiWaitbar(['Creating 4D Array y = ' num2str(j)],i/fL);
        end
        %waitbar(i/(fL),b,'Creating 4D array');
        % fprintf('.');
        if j < sL
            multiWaitbar(['Creating 4D Array y = ' num2str(j)],'Relabel',['Creating 4D Array y = ' num2str(j+1)]);
        end
    end
else
    
    for j = 1:sL
        for i = 1:fL
            if new
                if ind
                    HF1(i,j,:,:,:) = Read_Data2(loc,(i)+(fL*(j-1)),param,a);
                else
                    HF1(i,j,:,:) = Read_Data2(loc,(i)+(fL*(j-1)),param,a);
                end
            else
                [~,HF1(i,j,:,:)] = read_ucsdi_data(loc,(i)+(fL*(j-1)));
            end
            multiWaitbar(['Creating 4D Array y = ' num2str(j)],i/fL);
        end
        if j<sL
            multiWaitbar(['Creating 4D Array y = ' num2str(j)],'Relabel',['Creating 4D Array y = ' num2str(j+1)]);
        end
    end
end


    
    
    
%     for j = 1:sL
%         for i = 1:fL
%             if new
%                 [HF1{i,j}] = Read_Data2(loc,(i)+(fL*(j-1)),param); % Gets data sequentially
%             else
%                 [~,HF1{i,j}] = read_ucsdi_data(loc,(i)+(fL*(j-1))); % Gets data sequentially
%             end
%             % disp((i)+(fL*(j-1)));
%             multiWaitbar(['Creating 4D Array y = ' num2str(j)],i/fL);
%         end
%         %waitbar(i/(fL),b,'Creating 4D array');
%         % fprintf('.');
%         if j < sL
%             multiWaitbar(['Creating 4D Array y = ' num2str(j)],'Relabel',['Creating 4D Array y = ' num2str(j+1)]);
%         end
%     end
% end
% for i = 1:size(HF1,1)
%     for j = 1:size(HF1,2)
%     y(i,:,:) = HF1{i,1};
%     y2(j,:,:) = HF1{75,j};
%     end
% end

  

% if length(size(HF1{1})) > 2
%     if exist('a','var')
%         for i = 1:size(HF1,1)
%             for j = 1:size(HF1,2)
%                 HF1{i,j} = HF1{i,j}(:,:,a);
%             end
%         end
%     end
% end

% HF=zeros(size(HF1,1),size(HF1,2),size(HF1{1},1),size(HF1{1},2));
% for i = 1:fL
%     for j = 1:sL
%        % HF((i-1)*sL+j,:,:) = HF1{i,j}; %converts cells to pseudo 4-D array
%        HF(i,j,:,:) = HF1{i,j}; %Converts cell array to double
%     end
%     waitbar(.5+i/(fL*2),b,'Finalizing 4D array construction');
% end

%HF = reshape(HF,fL,sL,param.daq.HFdaq.pts,param.daq.HFdaq.NoBurstTriggers); %converts pseudo 4-D to 4-D
HF = 0;
if param.velmex.FastAxis == 'Y'
    HF = permute(HF,[2 1 3 4]); %rearrange Y and X if data was taken that way
end
%delete(b);
end



