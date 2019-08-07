%Input voltage of hydrophone signal in mV
%Output signal is in MPa
function P = mv2pressure(mv,b)

if exist('b','var')
    b = b;
else 
    b = 1;
end
if b
    P = mv/3.4087;
else
    P = mv/3.4087;
end