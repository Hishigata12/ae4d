function [slowpt, fast, Xstep, Ystep] = get_points(floc)
param = read_ucsdi_info(floc);
Tpts = param.velmex.XNStep*param.velmex.YNStep;
velmex = param.velmex;
slowpt = 1;
if velmex.FastAxis == 'X'
%     for i = 1:Tpts/XNStep
%     slowpt(i) = 
    slowpt = [slowpt (1:(velmex.YNStep-1))*velmex.XNStep];
else
    slowpt = [slowpt (1:(velmex.XNStep-1))*velmex.YNStep];
end
fast = velmex.FastAxis;
Xstep = param.velmex.XNStep;
Ystep = param.velmex.YNStep;

%param = read_ucsdi_info(floc);
