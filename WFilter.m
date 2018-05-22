% Creates a weiner filter for reducing noise
% Generally not useful, but might be helpful for removing noise in pulse
% data
% R asks for time axis type
% ln is line of M-Mode data to acquire noise profile

function [X] = WFilter(LF,TX)

if R == 'ae' || R == 'pe' || R == 'f' || R == 'fast' || R == 'M'