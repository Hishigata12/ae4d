%Input voltage of hydrophone signal in mV
%Output signal is in MPa
function P = mv2pressure(mv,b)

if exist('b','var')
    b = b;
else 
    b = 1;
end
if b
    P = 0.866*mv;
else
    P = 0.866*mv;
end