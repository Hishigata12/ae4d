function [pe, ae, lf] = signal_data()
global a
global f_loc
global xy
[a, HF, LF] = read_ucsdi_data(f_loc,xy);
pe{xy} = HF(:,:,1);
ae{xy} = HF(:,:,2);
lf{xy} = LF;
end
