%Outputs the LF and HF data of an AE scan
%Inputs: fname is the name of the info.txt file
%        n is the scan point number to grab

function [HF, LF] = Read_Data(fname,ScanPt,param,p)

if ~exist('ScanPt')
    ScanPt = 1;
end
ind = param.post.ind;
sweep = param.sweep;
if ind == 1 || sweep == 1
    froot = fname(1:end-9);
    fparam = [froot '_info.mat'];
    fLF = [froot '_LF.dat'];
    fHF = [froot '_HF.dat'];
else
    froot = fname(1:end-9);
    fparam = [froot '_info.mat'];
fLF = [froot '_LF_Avg.dat'];
fHF = [froot '_HF_Avg.dat'];
end

% temp = load(fparam);
% param = temp.bScanParm;

if ScanPt ~= 0
    
    % Gets LF data
    fid = fopen(fLF,'rb');
    nBytes = fread(fid,1,'int32');
    nPos = ftell(fid);
    stype = strtrim(fgets(fid,nBytes)); % LVDATA
    
    nver = fread(fid,1,'int32');
    LF = read_lfdata_v161108(fid,ScanPt,param);
    fclose(fid);
    
    
    %%%%%%%%%%%
    
    % Get HF Data
    if ~param.lf_only
        m = ScanPt;
        HF = read_hfdata(fHF,m,param,p);
        x = 1;
    else
        HF = [];
    end
    
    %%%%%%%%%%%
    
end
end

%%% FURTHER FUNCTIONS FOR LF AND HF %%%
function [data,param] = read_lfdata_v161108(fid,ScanPt,param)
n = fread(fid,1,'int32');
dsize = fread(fid,[1,n],'int32');

cur_pos = ftell(fid);

% first trace
param.nPts = dsize(2);
param.numChan = dsize(1);

% blk_size = (length(dsize) + 1)*4 + prod(dsize)*4;
blk_size = prod(dsize)*4;
offset_n = cur_pos + blk_size * (ScanPt-1);
fseek(fid,offset_n,'bof');


if param.post.ind
    data2 = fread(fid,fliplr(dsize),'single');
    dec = param.Scan.Avg;
    leng = size(data2,1)/dec;
    numchans = dsize(1);
    data2 = reshape(data2,[prod(dsize),1]);
    for i = 1:dec*numchans
        SingleAvg(:,i) = data2(leng*(i-1)+1:leng*i);
    end
    for i = 1:numchans
        for j = 1:dec
            data(:,j,i) = SingleAvg(:,i+(j-1)*(numchans));
        end
    end
elseif param.sweep
        data2 = fread(fid,fliplr(dsize),'single');
    dec = param.Scan.Avg;
    leng = size(data2,1)/dec;
    numchans = dsize(1);
    data3 = reshape(data2,[prod(dsize),1]);
    for i = 1:dec*numchans
        SingleAvg(:,i) = data2(leng*(i-1)+1:leng*i);
    end
    numchans = 2; %Temp Hardcode
    dec = dec/2; %Temp Hardcode
    for i = 1:numchans
        for j = 1:dec
            data(:,j,i) = SingleAvg(:,i+(j-1)*(numchans));
        end
    end
    data = squeeze(mean(data,2));
else
    data = fread(fid,fliplr(dsize),'single');
end
end


function [hfdata,param] = read_hfdata(hffile,ScanPt,param,p)
blk_idx = ScanPt;

if exist(hffile,'file')
    fid = fopen(hffile,'rb');
    blk_idx = 1;
else
    hf_avg_file = regexprep(hffile,'_P[0-9]{4,4}_','_');
    fid = fopen(hf_avg_file,'rb');
end

if fid < 0
    error('Cannot open HF data file');
end

sztype = fgets(fid,fread(fid,1,'int32'));
if(strcmpi(sztype(1:5),'ucsdi') == 0)
    fclose(fid);
    return;
end

nver = fread(fid,1,'int32');
switch(nver)
    case {140908,160714}
        [hfdata,param] = read_hfdata_v140908(fid);
    case {161108,170103}
        [hfdata,param] = read_hfdata_v161108(fid,ScanPt,param,p);
    case 170427
        [hfdata,param] = read_hfdata_v170427(fid,ScanPt);
    otherwise
        [hfdata,param] = read_hfdata_default(fid);
end
param.Version = nver;
fclose(fid);
end

function [data,parm] = read_hfdata_v161108(fid,ScanPt,param,p)

% nPos = ftell(fid);
n = fread(fid,1,'int32');
dsize = fread(fid,[1,n],'int32');

% blk_size = (length(dsize) + 1)*4 + prod(dsize)*4;
blk_size = prod(dsize)*4;
fseek(fid,blk_size*(ScanPt-1),'cof');
data = fread(fid,[prod(dsize),1],'single');
if param.post.ind || param.sweep
    tsize = dsize(1)*dsize(2)*dsize(3)*dsize(4); %Slow Time, Channels, Fast Time, Averages
else
    tsize = dsize(1)*dsize(2)*dsize(3);
end
if tsize ~= length(data)
    data = zeros(tsize,1);
    disp(['Missing point ' num2str(ScanPt)]);
end
if param.post.ind
    data = reshape(data,[dsize(3),dsize(2),dsize(1),dsize(4)]);
    data = permute(data,[3,1,4,2]);
    data = data(:,:,:,p);
elseif param.sweep
    d2 = reshape(data,floor(param.Daq.HF.Samples),param.Scan.Sweep_Tpt*param.Scan.Xpt,param.Scan.Avg);
    for i = 1:param.Scan.Xpt
        for j = 1:param.Scan.Sweep_Tpt
            HF(i,:,j,:) = d2(:,i+(param.Scan.Xpt*(j-1)),:);
        end
    end
    data = mean(HF,4);
    data = permute(data,[1 4 2 3]);
    
else
    data = reshape(data,fliplr(dsize));
    
    data = permute(data,[1,3,2]);
    data = data(:,:,p);
end

parm.dsize = dsize;
parm.ScanPt = ScanPt;
end
%
% figure;
% for i = 1:size(d2,3)
%     imagesc(squeeze(d2(:,:,i)))
%     pause(0.5);
% end
% for i = 1:param.Scan.Xpt
%     for j = 1:param.Scan.Sweep_Tpt
%         HF(i,:,j,:) = d2(:,i+(param.Scan.Xpt*(j-1)),:);
%     end
% end
% 
% figure; for i = 1:size(HF,3)
%     imagesc(HF2(:,:,i))
%     pause(0.4)
% end

