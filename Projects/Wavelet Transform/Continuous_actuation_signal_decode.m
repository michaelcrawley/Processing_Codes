function [ret_val,st_ids,F_f] = Continuous_actuation_signal_decode(sig_act,info_in)

%   Decodes the actuation times from continuous ramp signal that is a
%   surrogate of the actuation pulse train
%
%   The actuation signal for the LAFPAs is a rectangular pulse train (of
%   voltages).
%   Since the sampling rate of this signal is not high enough, the phase of
%   actuation cannot be unambiguously determined from sampling this raw
%   pulse train.
%   Instead, an arbitrary waveform generator is used to trigger a ramp
%   signal at each rising edge of the actuation pulse.
%   This ramp signal is high (at info_in.actuation_ramp.Vpp) by default,
%   but falls to 0 at the rising edge of the LAFPA control pulse train, and
%   then rises to the high value over a short time interval (which need not
%   be specified to the program). It then stays high until the next rising
%   edge of the actuation signal.
%   
%   Assumptions: The parameters of the ramp signal and the sampling are
%   chosen such that
%   1:  The slanted rising edge of each ramp is sampled at least thrice,
%       and the high value in between two ramps is sampled at least once.
%   2:  The falling edge of each ramp would be reliably detected by a
%       decrease in sampled voltages that is greater than Vpp/10.
%   3:  The actuation computer runs on a high-speed clock with rate
%       info_in.actuation_ramp.clk_rate. The actuation frequency is always
%       quantized such that the time period corresponds to an integer
%       number of clock pulses. The number of such clock pulses between
%       each ramp signal is determined here, and if they are not found to
%       be all the same, then an error is reported.
%
%   Input Parameters:
%   sig_act:	The sampled ramp signal. Dims are
%               1:  Block of data
%               2 to end:   Any other parameter (like different blocks)
%   info_in:    Struct that must contain the following metadata fields
%               at least
%       sampling_rate:  Sampling rate (same as rate of pressure sampling)
%       actuation_ramp: Struct with fields
%           Vpp:        High value of ramp signal (in Volts)
%           clk_rate:   Actuation controller base clock rate (in Hz)
%           period:     Rise time period
%
%   Output Parameters:
%   ret_val:Error code
%       0:  no error
%       -n: Unambiguous # clock samples couldn''t be detected in block #n
%       1:  All blocks of data did not have same forcing rate
%   st_ids:	Cell matrix of dimensions equal to the 2nd to last dimensions
%           of sig_act. Each cell contains the index numbers of the start
%           of the pulses in the block of data. The index numbers may not
%           be integers, since the beginnings of the pulses may not
%           coincide with the sampling instants.
%   F_f:    Actuation frequency as determined here (in Hz)

%   Correct processing indicator
ret_val = 0;

%   Threshold for determining when a falling edge of the RAMP signal is
%   encountered
threshold = info_in.actuation_ramp.Vpp/10;

%   Columnize sig_act if it is a vector
if isvector(sig_act)
    sig_act = sig_act(:);
end

%   Determine the block size (assumed to be the size of first dim)
sig_sz = size(sig_act);
blk_sz = sig_sz(1);
%   Pack all dims of sig_act beyond the 1st dim into the 2nd dim
sig_act = reshape(sig_act,blk_sz,[]);
%   Determine the number of separate blocks (each having a single PIV laser
%   fire)
n_blks = size(sig_act,2);

%   Increase in ramp voltage over a sampling period
rise_rate = info_in.actuation_ramp.Vpp/(info_in.actuation_ramp.period*info_in.sampling_rate);

%   Pre-allocate and initialize default return values
st_ids = cell(1,n_blks); 

N_Clk_samples_all = zeros(1,n_blks); %Pre-allocate
for ndx = 1:n_blks
    signal_n = sig_act(:,ndx);
    
    %   Determine the indices where the signal decreases significantly
    start_idx = find(diff(signal_n) < -threshold);

    %   The falling edges are sharp, so that at most one sample can
    %   coincide with a falling edge. 
    %   However, the immediate succeeding sample may be lower than such a
    %   sample on the falling edge, even though it is on the rising ramp.
    %   To pick off the rising edge unambiguosly, reject the first decrease
    %   if the immediate next behaviour is also a decrease
    start_idx = setdiff(start_idx,start_idx - 1);

    %   Owing to the behavior of the 'diff' function, the actual critical
    %   indices (which have the voltage lower than the 'high' value) are
    %   the successors of the ones identified above.
    %   However, at times, even such succeding samples may be on the
    %   falling edges of the ramps.
    %   Thus, using the guarantee that at least 2 samples will be on the
    %   ramps, the succeding samples are to be used as sure to be on the
    %   ramp.
    %   The following is then the array of indices that are guaranteed to
    %   be on the ramps
    ramp_idx = start_idx + 2;

    %   If the last 'ramp_idx' falls outside the range of block samples,
    %   then discard it
    if ramp_idx(end) > blk_sz%-1
        ramp_idx = ramp_idx(1:end-1);
    end

    %   The sample index that coincides with the zero phase of each ramp
    %   would not be an integer in general.
    %   This is addressed in the current implementation.
    %   First determine such (possibly non-integer) indices.
    %   For this, the voltage on the ramp, and the known voltage rise per
    %   sample on the ramp are used.
    st_ids{ndx} = ramp_idx(:) - sig_act(ramp_idx(:))/rise_rate;

    
%     ss_rise_1st = sig_act(ramp_idx,ndx); %Voltages at first point in rising ramp
%     ss_rise_2nd = sig_act(ramp_idx + 1,ndx); 
%     st_ids{ndx} = ramp_idx - ss_rise_1st./(ss_rise_2nd - ss_rise_1st);
    
    %   Determine the time interval between each actuator pulse
    T_act_all = diff(st_ids{ndx})/info_in.sampling_rate;
    
    %   Determine # clock samples
    N_Clk_samples = round(T_act_all*info_in.actuation_ramp.clk_rate);
    
    %   Ensure that all # clock samples come out as identical
    if norm(diff(N_Clk_samples))
        ret_val = -ndx;
    end
    
    %   Note the common N_clk_samples for this block
    N_Clk_samples_all(ndx) = N_Clk_samples(1);
end

%   Reshape st_ids to original size of sig_act
st_ids = reshape(st_ids,[1,sig_sz(2:end)]);

%   Ensure that # clock samples come out as identical for all blocks
if norm(diff(N_Clk_samples_all))
    ret_val = 1;
end

%   Determine unique actuation frequency
F_f = info_in.actuation_ramp.clk_rate/N_Clk_samples_all(1);
