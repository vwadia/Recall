function [ppc, f, ppc_boot] = compute_ppc(lfp_data,spike_data,min_freq,max_freq,n_freq,whichSteps,fs,run_boot)
% Function gifted to me by Jonathan
% December2022
% input:

% lfp_data: samples x channels X trials raw data (samples = indices based on sampling rate)
% spike_data: samples x neurons X trials raw data (same sampling rate and time as lfp_data - note samples are 0s and 1s like raster)
% low_freq: lowest frequency to compute
% high_freq: highest frequency to compute
% n_freq: how many frequency steps?
% whichSteps: 'linear' or 'log'
% fs: sampling rate
% run_boot - to compute bootstrap distribution (added by Varun)

% output:
% ppc: neurons x channels x freq
% f: computed center frequencies
if strcmp(run_boot, 'false')
    run_boot = 0;
elseif strcmp(run_boot, 'true')
    run_boot = 1;
end

n_trials = size(lfp_data,3);
n_channels = size(lfp_data,2);
n_neurons = size(spike_data,2);

% Wavelet convolution parameters
if strcmp(whichSteps,'log')
    f = logspace(log10(min_freq),log10(max_freq),n_freq);
elseif strcmp(whichSteps,'linear')
    f = linspace(min_freq,max_freq,n_freq);
end

wl_time      = -2:1/fs:2;
wl_time_half = (length(wl_time)-1)/2;

% FFT parameters (use next-power-of-2)
n_samples_wl          = length(wl_time);
n_samples_data        = size(lfp_data,1);
n_samples_convolution = n_samples_wl+n_samples_data-1;
n_samples_conv_pow2   = pow2(nextpow2(n_samples_convolution));
if strcmp(whichSteps,'log')
    wavelet_cycles = logspace(log10(3),log10(10),n_freq);
elseif strcmp(whichSteps,'linear')
    wavelet_cycles = linspace(3,10,n_freq);
end

% Computation and FFT of wavelets
% fprintf('Computation and fft of wavelets for all %d specified frequencies...\n',n_freq)
% compute wavelets for each frequency
wavelets = zeros(n_samples_wl,n_freq);
for fi = 1:n_freq
    wavelets(:,fi) = (pi*f(fi)*sqrt(pi))^-.5 * exp(2*1i*pi*f(fi).*wl_time) .* exp(-wl_time.^2./(2*( wavelet_cycles(fi) /(2*pi*f(fi)))^2))/f(fi); % first term is confusing, 
 end

% Fourier transform of the wavelets
wavelets_fft = fft(wavelets, n_samples_conv_pow2);

phase_data = zeros(n_samples_data, n_trials, n_channels, n_freq);

%% Wavelet convolution
lastsize = 0;
for i_trial = 1:n_trials
    % Tell which trial is processed
    fprintf(repmat('\b', 1, lastsize));
%     lastsize = fprintf('Processing trial %d of %d',i_trial,n_trials);
    
    % cut out this trial
    this_trial = squeeze(lfp_data(:,:,i_trial));
    
    % FFT of data (note: this doesn't change on frequency iteration)
    fft_trial = fft(this_trial,n_samples_conv_pow2);
    
    % compute convolution for each frequency
    for i_freq=1:n_freq
        
        % duplicate wavelets to match number of channel
        wl = repmat(wavelets_fft(:,i_freq), [1 n_channels]);
        
        % run convolution
        convResult = ifft(wl.*fft_trial,n_samples_conv_pow2);
        convResult = convResult(1:n_samples_convolution,:); % here the extra points from the power-of-2 FFT are removed
        convResult = convResult(wl_time_half+1:end-wl_time_half,:);
        
        % Put averaged data to tf-matrix
        phase_data(:,i_trial,:,i_freq) = angle(convResult);
        
    end %freq
end %trial
fprintf('\n')

%% Concatenate trials
phase_data = reshape(phase_data, [n_samples_data * n_trials n_channels n_freq]);

%% Compute PPC
n_iter = 500;
ppc_boot = nan(n_neurons,n_channels,n_freq,n_iter);
ppc = nan(n_neurons,n_channels,n_freq);
lastsize = 0;
for ineuron = 1:n_neurons
    fprintf(repmat('\b', 1, lastsize));
%     lastsize = fprintf('Computing PPC in neuron %d of %d',ineuron,n_neurons);
    
    this_spike_data_wTrial = squeeze(spike_data(:,ineuron,:));
    this_spike_data = logical(reshape(this_spike_data_wTrial, [n_samples_data * n_trials 1]));
    
    n_spikes = sum(this_spike_data);

    sfc_neuron  = squeeze(exp(1i*phase_data(this_spike_data,:,:)));
    outsum = squeeze(nansum(sfc_neuron));
    ppc(ineuron,:,:) = (outsum.*conj(outsum) - n_spikes)./(n_spikes*(n_spikes-1));

    % Compute surrogate distribution
    if run_boot 
        fprintf(repmat('\b', 1, lastsize));
%         lastsize = fprintf('Computing PPC distribution for neuron %d of %d',ineuron,n_neurons);
        for n_rep = 1:n_iter
            spike_data_boot = Utilities.Shuffle(this_spike_data_wTrial, 2); % shuffle trials
            spike_data_boot = logical(reshape(spike_data_boot, [n_samples_data * n_trials 1]));

            sfc_neuron_boot  = squeeze(exp(1i*phase_data(spike_data_boot,:,:)));
            outsum_boot = squeeze(nansum(sfc_neuron_boot));
            ppc_boot(ineuron,:,:,n_rep) = (outsum_boot.*conj(outsum_boot) - n_spikes)./(n_spikes*(n_spikes-1));
        end
    end
end
fprintf('\n')




