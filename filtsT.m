%X is filtered output
%x is input
%p are the window sizes or interp factor
%use declares whether a filtering method will be considered
%param holds parameters (not used)

function X = filtsT(x,p,use,param)
if size(x) == 2
    x = permute(x,[1 4 3 2]);
elseif size(x) == 3
    x = permute(x,[1,4,2,3]);
end

if use(1)
    if p(1) > 1
        if mod(p(1),2) == 0
            errordlg('window must be odd')
        else
            for i = 1:size(x,1)
                for j = 1:size(x,2)
                    for k = 1:size(x,3)
                        x2(i,j,k,:) = medfilt1(squeeze(x(i,j,k,:)),p(1));
                    end
                end
                multiWaitbar('Median Filtering in Time',i/size(x,1));
            end
        end
    end
    x = x2;
    clear x2
end

if use(2)
    if p(2) > 1
        if mod(p(2),2) == 0
            errordlg('window must be odd')
        else
            a = 1;
            b = (1/p(2))*ones(1,p(2));
            for i = 1:size(x,1)
                for j = 1:size(x,2)
                    for k = 1:size(x,3)
                        x2(i,j,k,:) = filter(b,a,squeeze(x(i,j,k,:)));
                    end
                end
                multiWaitbar('Smoothing in Time',i/size(x,1));
            end
        end
    end
    x = x2;
    clear x2
end

if use(3)
    if p(3) < 1
        Q = 'Decimating';
    else
        Q = 'Interpolating';
    end
    if p(3) ~= 1
        for i = 1:size(x,1)
            for j = 1:size(x,2)
                for k = 1:size(x,3)
                    x2(i,j,k,:) = interp1(linspace(0,1,size(x,4)),squeeze(x(i,j,k,:)),linspace(0,1,size(x,4)*p(3)));
                end
            end
            multiWaitbar([Q ' in Time'],i/size(x,1));
        end
        
    end
    x = x2;
    clear x2
end

multiWaitbar('CLOSEALL');

X = x;