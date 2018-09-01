%Outputs the LF and HF data of an AE scan
%Inputs: fname is the name of the info.txt file
%        n is the scan point number to grab

function [LF, HF] = Read_Data(fname,ScanPt)

if ~exist('ScanPt')
    ScanPt = 1;
end

froot = fname(1:end-4);
fparam = [froot '_info.mat'];
fLF = [froot '_LF_Avg.dat'];
fHF = [froot '_HF_Avg.dat'];

temp = load(fparam);
param = temp.bScanParm;

if ScanPt ~= 0
    
    % Gets LF data
    fid = fopen(fLF,'rb');
    nBytes = fread(fid,1,'int32');
    nPos = ftell(fid);
    stype = strtrim(fgets(fid,nBytes)); % LVDATA
    
    nver = fread(fid,1,'int32');
    LF = read_lfdata_v161108(fid,ScanPt);
    
    
    %%%%%%%%%%%
    
    % Get HF Data
    HF = read_hfdata(fHF,ScanPt);
    
    %%%%%%%%%%%
    
end
end

%%% FURTHER FUNCTIONS FOR LF AND HF %%%
function [data,param] = read_lfdata_v161108(fid,ScanPt)
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
end


function [hfdata,param] = read_hfdata(hffile,ScanPt)
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
        [hfdata,param] = read_hfdata_v161108(fid,blk_idx);
    case 170427
        [hfdata,param] = read_hfdata_v170427(fid,ScanPt);
    otherwise
        [hfdata,param] = read_hfdata_default(fid);
end
param.Version = nver;
fclose(fid);
end

function [data,parm] = read_hfdata_v161108(fid,ScanPt)

% nPos = ftell(fid);
n = fread(fid,1,'int32');
dsize = fread(fid,[1,n],'int32');

% blk_size = (length(dsize) + 1)*4 + prod(dsize)*4;
blk_size = prod(dsize)*4;
fseek(fid,blk_size*(ScanPt-1),'cof');
data = fread(fid,[prod(dsize),1],'single');
tsize = dsize(1)*dsize(2)*dsize(3);
if tsize ~= length(data)
    data = zeros(tsize,1);
    disp(['Missing point ' num2str(ScanPt)]);
end
data = reshape(data,fliplr(dsize));

data = permute(data,[1,3,2]);

parm.dsize = dsize;
parm.ScanPt = ScanPt;
end


