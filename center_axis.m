function center_axis(Sig,Ax,Axis,a,scaling)
%Axis must be 'HF' or 'LF'
%a must be length 2 vector for both min and max of axis
%use 'lin' for linear scale, otherwise it will use dB

%global AE_axis;
%global m_axis;
%global AE_f;
%global AE_nf;

if Axis == 'HF'
    %ax = 'ylim';
    a2(1) = find(Ax>a(1),1);
    a2(2) = find(Ax>a(2),1);
elseif Axis == 'LF'
    %ax = 'xlim';
    a2(1) = find(Ax>a(1),1);
    a2(2) = find(Ax>a(2),2);
else
    errordlg('Axis must be HF or LF','You are ass hole');
end

if Axis == 'HF'
    ylim(a);
elseif Axis =='LF'
    xlim(a);
end

if scaling  == 'lin' 
axmax = max(Sig(a2(1):a2(2)));
axmin = min(Sig(a2(1):a2(2)));
else
end
caxis([axmin axmax])

end
        
