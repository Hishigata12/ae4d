function [daq velmex hp] = testcode(b);
%This must exist at the beginning of all codes
global floc
[param HF LF] = read_ucsdi_data(floc,b);
velmex = param.velmex;
hp = param.hp;
daq = param.daq;