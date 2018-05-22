%Get parameters
[f p]  = uigetfile(fullfile(pwd,'*_info.dat'));
[param] = read_ucsdi_info([p f]);

%Get File name
[file path] = uigetfile(fullfile(pwd,'*.bsq'));
fn = [path file];

%Get size matrix
ro = param.daq.HFdaq.pts;
co = param.velmex.XNStep;
ba = param.velmex.YNStep;
sz = [ro, co, ba];

%precision
prec = 'double';

%declare offset for where data begins
offs = 100;

%interleave
il = 'bsq';

%endian
en = 'ieee-be';

a = multibandread(fn,sz,prec,offs,il,en);