function [M, ax] = make_axes(param,dims,delay,D,t)
if ~exist('delay','var')
    delay = 5.2;
end
M.Depth = 1.48*param.daq.HFdaq.pts/param.daq.HFdaq.fs_MHz-delay;  % Max Depth
M.x = linspace(0,M.Depth,param.daq.HFdaq.pts);
M.stime = param.daq.HFdaq.duration_ms; %duration in of slow time
ax.depth = linspace(-delay,M.Depth,dims(3));
ax.stime = linspace(0,M.stime,dims(4));
if param.velmex.XDist ~= 0
    if param.velmex.SlowAxis == 'X'
         ax.x = linspace(-param.velmex.YDist/2,param.velmex.YDist/2,dims(1));
    else
    ax.x = linspace(-param.velmex.XDist/2,param.velmex.XDist/2,dims(1));
    end
else
    ax.x = 1;
end
if param.velmex.YDist ~= 0
     if param.velmex.SlowAxis == 'X'
          ax.y = linspace(-param.velmex.XDist/2,param.velmex.XDist/2,dims(2));
     else
    ax.y = linspace(-param.velmex.YDist/2,param.velmex.YDist/2,dims(2));
     end
else
    ax.y = 1;
end

if exist('D','var')
    if exist('t','var')
        M.xL = find(ax.depth < D(2));
        M.xH = find(ax.depth > D(1));
        M.xT = intersect(M.xL,M.xH);
        M.sT = find(ax.stime > t,1);
    end
end
end