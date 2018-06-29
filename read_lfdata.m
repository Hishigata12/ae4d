function [LFData,param] = read_lfdata(f_root,ScanPt)

parascan = evalin('base','parascan');
f_path = parascan.f_path;

NX = parascan.velmex.XNStep;
NY = parascan.velmex.YNStep;

nScanPts = NX * NY;
if exist('ScanPt','var') == 0
    ScanPt = 1:nScanPts;
end

numRead = length(ScanPt);
for ni = 1:numRead
    fname = fullfile(f_path,f_root);
    [data,param] = read_lfdata_1(fname,ScanPt(ni));
    if ni == 1 && numRead > 1
        LFData = repmat(data,[ones(1,ndims(data)),numRead]);
    end
    LFData(:,:,ni) = data;
end

%%%%%%%%%%%%%%%%%%%%%%%%
% 
% if nver >= 161108
%     [LFData] = read_lfdata(fullfile(f_path,lffile),ScanPt);
% else
%     
%     fid = fopen(fullfile(f_path,lffile),'rb');
%     if(fid < 0)
%         uiwait(errordlg([lffile,' does not exist'],mfilename,'modal'));
%         return;
%     end
%     
%     try
%         n = fread(fid,1,'int32');
%         lfsize = fread(fid,[1,n],'int32');
%         LFData = fread(fid,fliplr(lfsize),'single');
%         fclose(fid);
%     catch ME
%         fclose(fid);
%         uiwait(errordlg(ME.message,mfilename,'modal'));
%         return;
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%


end

function [LFData,param] = read_lfdata_1(lffile,ScanPt)

if exist([lffile,'_LF_Avg.dat'],'file')
    fid = fopen([lffile,'_LF_Avg.dat'],'rb');
    blk_idx = ScanPt;
else
    lf_file = [lffile,sprintf('_P%04d_LF_Avg.dat',ScanPt)];
    fid = fopen(lf_file,'rb');
    blk_idx = 1;
end

if fid < 0
    error('Cannot open LF data file');    
end

szFileType = 'ucsdi_lf';

nBytes = fread(fid,1,'int32');
nPos = ftell(fid);
if(nBytes < 1 || nBytes > 64)
    fclose(fid);
    error('Invalid data file');
end

stype = strtrim(fgets(fid,nBytes)); % LVDATA

if(strcmpi(stype,szFileType)==1)
    nver = fread(fid,1,'int32');
    switch nver
        case {161108,170103}
            [LFData,param] = read_lfdata_v161108(fid,blk_idx);
        case 170427
            [LFData,param] = read_lfdata_v170427(fid,blk_idx);
        case 170612
            [LFData,param] = read_lfdata_v170612(fid,blk_idx);
        otherwise
            LFData = [];
            param = struct;
    end
else
    fseek(fid,nPos,'bof');
    lfsize = fread(fid,[1,nBytes],'int32');
    LFData = fread(fid,fliplr(lfsize),'single');
%     fclose(fid);
%     LFData = [];
    param.nPts = lfsize(2);
    param.numChan = lfsize(1);
end

fclose(fid);

end


function [data,param] = read_lfdata_v161108(fid,ScanPt)    
% npos = ftell(fid);

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

data = fread(fid,fliplr(dsize),'single');

% npos1 = ftell(fid);
% fseek(fid,0,'eof');
% nlen = ftell(fid);
% 
% numTrace = floor((nlen - npos)/(npos1 - npos));
% fseek(fid,npos1,'bof');
% 
% data = repmat(data,[ones(1,length(dsize)),numTrace]);
% for ni=2:numTrace
%     dsize = fread(fid,[1,fread(fid,1,'int32')],'int32');
%     data(:,:,ni) = fread(fid,dsize,'single');
% end
% 
% param.NumAcq = numTrace;

end


function [data,param] = read_lfdata_v170427(fid,ScanPt)    
% npos = ftell(fid);

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

data = fread(fid,fliplr(dsize),'single');

% npos1 = ftell(fid);
% fseek(fid,0,'eof');
% nlen = ftell(fid);
% 
% numTrace = floor((nlen - npos)/(npos1 - npos));
% fseek(fid,npos1,'bof');
% 
% data = repmat(data,[ones(1,length(dsize)),numTrace]);
% for ni=2:numTrace
%     dsize = fread(fid,[1,fread(fid,1,'int32')],'int32');
%     data(:,:,ni) = fread(fid,dsize,'single');
% end
% 
% param.NumAcq = numTrace;

end


function [data,param] = read_lfdata_v170612(fid,ScanPt)    
% npos = ftell(fid);

n = fread(fid,1,'int32');
dsize = fread(fid,[1,n],'int32');

cur_pos = ftell(fid);

% first trace
param.nPts = dsize(end);
param.numChan = dsize(end-1);

% blk_size = (length(dsize) + 1)*4 + prod(dsize)*4;
blk_size = prod(dsize)*4;
offset_n = cur_pos + blk_size * (ScanPt-1);
fseek(fid,offset_n,'bof');

data = fread(fid,fliplr(dsize),'single');

% npos1 = ftell(fid);
% fseek(fid,0,'eof');
% nlen = ftell(fid);
% 
% numTrace = floor((nlen - npos)/(npos1 - npos));
% fseek(fid,npos1,'bof');
% 
% data = repmat(data,[ones(1,length(dsize)),numTrace]);
% for ni=2:numTrace
%     dsize = fread(fid,[1,fread(fid,1,'int32')],'int32');
%     data(:,:,ni) = fread(fid,dsize,'single');
% end
% 
% param.NumAcq = numTrace;

end
