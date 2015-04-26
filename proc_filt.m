function dat= proc_filt(dat, b, a)
%PROC_FILT - Digital filtering
%
%Usage:
% DAT= proc_filt(DAT, B, A)
%
% Apply digital (FIR or IIR) filter. This filtering is causal, but
% produces a phase shift.
%
%Input:
% DAT   - data structure of continuous or epoched data
% B,A   - filter coefficients
%
%Output:
% DAT   - updated data structure
%
%Example:
% % Let cnt be a structure of multi-variate time series ('.x', time along first
% % dimension) with sampling rate specified in field '.fs'.
% [b,a]= butter(5, [7 13]/cnt.fs*2);
% % Apply a causal band-pass filter 7 to 13Hz to cnt:
% cnt_flt= proc_filt(cnt, b, a);
%
%See also proc_filtfilt.

dat.x(:,:)= filter(b, a, dat.x(:,:));



% Dec 2008: copied from IDA toolbox. 
% All rights belong to the authors and Fraunhofer FIRST.IDA.
% http://ida.first.fraunhofer.de/homepages/ida/


