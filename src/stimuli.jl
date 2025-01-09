export pure_tone, gaussian_noise, gaussian_noise_sl, sam_noise, sam_noise_sl, te_noise

"""
    pure_tone(freq, phase, dur, fs)

Synthesize pure tone; `freq` (Hz), phase `ϕ` (rads), `dur` (s), sample rate `fs` (Hz)
"""
pure_tone(freq, ϕ, dur, fs) = sin.(2π .* freq .* timevec(dur, fs) .+ ϕ)
pure_tone(; freq=1e3, ϕ=0.0, dur=1.0, fs=100e3) = pure_tone(freq, ϕ, dur, fs)

"""
    gaussian_noise(; kwargs...)

Synthesize bandpass Gaussian noise at a specified overall level

Synthesizes bandpass Gaussian noise by first selecting a random vector of numbers from the
unit normal distribution and then filtering them with (by default) a 4th order Butterworth 
bandpass filter from `freq_low` to `freq_high` (although any filter can be passed in its
place via the `filter` argument, in which case `freq_low` and `freq_high` are ignored).

# Arguments:
- `freq_low=0.5e3`: Lower cutoff of the bandpass filter (Hz)
- `freq_high=20e3`: Upper cutoff of the bandpass filter (Hz)
- `level=50.0`: Overall level before ramping (dB SPL)
- `dur=1.0`: Total duration, including ramps (s)
- `dur_ramp=0.01`: Ramp duration (s)
- `fs=100e3`: Sampling rate (Hz)
- `filter`: Filter used to bandpass filter noise (by default, 4th order Butterworth)
"""
function gaussian_noise(;
    freq_low=0.05e3,
    freq_high=20e3,
    level=50.0,
    dur=1.0,
    dur_ramp=0.01,
    fs=100e3,
    filter=digitalfilter(Bandpass(freq_low, freq_high; fs=fs), Butterworth(4)),
)
    # Generate random numbers and filter
    waveform = filt(filter, randn(samples(dur, fs)))

    # Scale to overall level
    scale_dbspl!(waveform, level)

    # Ramp and return
    cosine_ramp(waveform, dur_ramp, fs)
end

"""
    gaussian_noise_sl(; kwargs...)

Same as `gaussian_noise`, except that `level` is interpreted as spectrum level
"""
function gaussian_noise_sl(;
    freq_low=0.05e3,
    freq_high=20e3,
    level=10.0,
    dur=1.0,
    dur_ramp=0.01,
    fs=100e3,
    filter=digitalfilter(Bandpass(freq_low, freq_high; fs=fs), Butterworth(4)),
)
    # Generate random numbers and scale to correct overall level based on spectrum level
    waveform = scale_dbspl(randn(samples(dur, fs)), sl_to_ol(level, fs))

    # Apply filter to waveform
    waveform = filt(filter, waveform)

    # Ramp and return
    cosine_ramp(waveform, dur_ramp, fs)
end

"""
    te_noise(; kwargs...)

Synthesize threshold-equalizing noise (Moore et al., 2000)

Synthesize threshold-equalizing noise (TEN), based on code adapted from from original code
from lab of A. J. Oxenham. TEN is originally described in:

Moore, B. C. J., Huss, M., Vickers, D. A., Glasberg, B. R., and Alcántra, J. I. (2000). “A
test for the diagnosis of dead regions in the cochlea,” British Journal of Audiology, 34,
205–224. doi:10.3109/03005364000000131

TEN is a noise spectrally shaped so as to acheive uniform masked thresholds for a pure tone
as a function of frequency from 125–15000 Hz. TEN is produced under the assumption that the
power of a signal at masked threshold is:

    P₀ = N₀ × K × ERB

where N₀ is the noise power spectral density, K is the signal-to-noise ratio at the output 
of the auditory filter required for masked threshold at the given frequency (this varies
with frequency), and ERB is the equivalent rectangular bandwidth of the auditory filter
at the given frequency.

Note that TEN level usually refers to level in the ERB centered at 1 kHz, and NOT overall
level. 

# Arguments:
- `freq_low=0.25e3`: Lower cutoff of the noise passband (Hz)
- `freq_high=15e3`: Upper cutoff of the noise passband (Hz)
- `level=20.0`: Level in the 1 kHz ERB before the ramp is applied (dB SPL)
- `dur=1.0`: Total duration, including ramps (s)
- `dur_ramp=0.01`: Ramp duration (s)
- `fs=100e3`: Sampling rate (Hz)
"""
function te_noise(;
    freq_low=0.05e3,
    freq_high=20e3,
    level=50.0,  # in the 1-kHz ERB!!
    dur=1.0,
    dur_ramp=0.01,
    fs=100e3,
)
    # Determine length of stimulus in samples
    n = samples(dur, fs)

    # Calculate some values we'll need. These values include `idxs`, which is a boolean
    # indexing vector indicating which of the frequency bins are included in the stimulus
    # passband, considering only positive frequencies below the Nyquist frequency.
    binsize = fs/n  # freq-bin resolution in Hz
    bin_low = Int(round(freq_low/binsize) + 1)
    bin_high = Int(round(freq_high/binsize) + 1)
    n_bin = Int(bin_high - bin_low + 1)
    idxs = bin_low:bin_high  # freq-bin indices included in the noise passband
    all_freqs = (0:(n-1)) .* binsize ./ 1e3  # freq-bin center freqs in kHz

    # Allocate empty arrays for real and imaginary coefficients of noise in spectral domain
    a = zeros(n)
    b = zeros(n)

    # Fill space between low and high bins with random N(0, 1) values. We will use a and b
    # as random real and imaginary coefficients to synthesize a noise in the spectral domain
    # (see next step)
    a[idxs] .= rand(Normal(0, 1), n_bin)
    b[idxs] .= rand(Normal(0, 1), n_bin)

    # Create stimulus spectrum by using a and b as real and complex Fourier coefficients
    fspec = a .+ 1im .* b

    # Set values of K 
    # These were obtained by personal correspondence from B. C. J. Moore to A. J. Oxenham 
    # via the reference source code; also reported in:
    #
    # Moore, B. C., Glasberg, B. R., & Baer, T. (1997). A model for the prediction of
    # thresholds, loudness, and partial loudness. Journal of the Audio Engineering Society,
    # 45(4), 224-240.
    #
    # Note that in the equation above, K is a linear ratio; here, K is recorded as a 
    # decibel quantity instead.
    K = [[0.0500,  13.5000],
         [0.0630,  10.0000],
         [0.0800,   7.2000],
         [0.1000,   4.9000],
         [0.1250,   3.1000],
         [0.1600,   1.6000],
         [0.2000,   0.4000],
         [0.2500,  -0.4000],
         [0.3150,  -1.2000],
         [0.4000,  -1.8500],
         [0.5000,  -2.4000],
         [0.6300,  -2.7000],
         [0.7500,  -2.8500],
         [0.8000,  -2.9000],
         [1.0000,  -3.0000],
         [1.1000,  -3.0000],
         [2.0000,  -3.0000],
         [4.0000,  -3.0000],
         [8.0000,  -3.0000],
         [10.0000,  -3.0000],
         [15.0000,  -3.0000]]
    K = hcat(K...)
    K_freq = K[1, :]  # the frequencies at which K is defined in kHz
    K = K[2, :]       # values of K at corresponding frequencies (dB)

    # Interpolate K over all_freqs using a cubic spline
    # NOTE: This step is surprisingly tricky and poorly documented/tested in other 
    # implementations of TEN. Here, we use a cubic spline interpolation where extrapolated
    # values are simply set to the endpoints (i.e., values for frequencies less than 0.05 
    # kHz are set to 13.5, values for frequencies greater than 15 kHz are set to -3.0). The
    # interpolation is performed on a log(kHz) axis.
    spline = CubicSpline(log.(K_freq), K; extrapl=0.0, extrapr=0.0)
    K_interp = spline[log.(all_freqs)]
    K_interp[1] = 13.5;  # set 0 Hz bin k_interp to match other low-freq bins, avoid NaN

    # Calculate ERB at each frequency bin
    ERB = 24.7 .* ((4.37 .* all_freqs) .+ 1)  # see: Eq. 3, Glasberg and Moore (1990), Hear. Res.

    # Calculate critical ratio, in dB overall level, at every frequency bin
    # Under the assumptions of this stimulus/power spectrum model of masking, the value
    # `cr_ERB` in each frequency bin is the overall level of the signal that would be 
    # just detectable at the corresponding frequency, assuming the noise were spectrally
    # flat and at 0 dB SPL spectrum level. This is explained helpfully by Patterson and 
    # colleagues in Section IV of:
    #
    # Patterson, R. D., Nimmo‐Smith, I., Weber, D. L., & Milroy, R. (1982). The
    # deterioration of hearing with age: Frequency selectivity, the critical ratio, the
    # audiogram, and speech threshold. The Journal of the Acoustical Society of America,
    # 72(6), 1788-1803.
    cr_ERB = K_interp .+ (10 .* log10.(ERB))
    TEN_No = -cr_ERB

    # Next, we spectrally shape the noise according to the critical ratio calculated above;
    # shaping is done on a bin-by-bin basis.
    fspec[idxs] .= fspec[idxs] .* 10 .^ (TEN_No[idxs] ./ 20)

    # Next, we calculate the level in the ERB centered at 1 kHz in the spectral domain;
    # this level estimate is used to scale the noise below
    index1kERB = (all_freqs .> 0.935) .& (all_freqs .< 1.0681)
    rms_1kHz = sqrt(1/n^2 * real(sum(fspec[index1kERB] .* conj(fspec[index1kERB])))) / sqrt(2)
    level_1kHz = 20*log10(rms_1kHz/20e-6)

    # Use IFFT to get the time-domain representation of the noise
    # Note that since we're implicitly working with the one-sided spectrum, we take the real
    # part to appropriately cast half the energy into the other side of the spectrum
    noise = real(ifft(fspec))

    # Scale noise; we amplify/attenuate the noise by the difference between the requested 
    # level (specified in terms of the level in the ERB centered at 1 kHz) and the empirical
    # level in the 1-kHz ERB
    noise = amplify(noise, level - level_1kHz)

    # Apply raised-cosine ramp as final step
    noise = cosine_ramp(noise, dur_ramp, fs)

    return noise

end


"""
    sam_noise(; kwargs...)

Synthesize sinusoidally amplitude-modulated noise with mod freq `fm` and depth `m`

Synthesizes a bandlimited Gaussian noise carrier using `gaussian_noise` (additional 
unspecified keyword arguments are passed to `gaussian_noise`, and can thus be used to 
control the bandwidth or filters used in synthesizing the noise) and then modulates it
using a sinusoidal modulator. Here, `level` is the overall level in dB SPL after modulation.

# Arguments
- `fm=10.0`: Modulation rate (Hz)
- `ϕ=-π/2`: Modulation phase (radians)
- `m=1.0`: Modulation depth (a.u.)
- `dur=1.0`: Total duration, including ramps (s)
- `dur_ramp=0.01`: Ramp duration (s)
- `level=50.0`: Overall level after modulation (dB SPL)
- `fs=100e3`: Sampling rate (Hz)
"""
function sam_noise(;
    fm=10.0,
    ϕ=-π/2,
    m=1.0,
    dur=1.0,
    dur_ramp=0.1,
    level=50.0,
    fs=100e3,
    kwargs...
)
    # Synthesize noise carrier
    carrier = gaussian_noise(; dur=dur, dur_ramp=dur_ramp, level=level, fs=fs, kwargs...)
    level_pre_mod = dbspl(carrier)

    # Synthesize modulator
    if fm == 0.0
        modulator = zeros(length(carrier))
    else
        modulator = pure_tone(fm, ϕ, dur, fs)
    end

    # Combine carrier and modulator and scale
    stim = (1.0 .+ m .* modulator) .* carrier
    scale_dbspl(stim, level_pre_mod)
end

"""
    sam_noise_sl(; kwargs...)

Same as `sam_noise`, except that `level` is interpreted as spectrum level 
"""
function sam_noise_sl(;
    fm=10.0,
    ϕ=-π/2,
    m=1.0,
    dur=1.0,
    dur_ramp=0.1,
    level=50.0,
    fs=100e3,
    kwargs...
)
    # Synthesize noise carrier
    carrier = gaussian_noise_sl(; dur=dur, dur_ramp=dur_ramp, level=level, fs=fs, kwargs...)
    level_pre_mod = dbspl(carrier)

    # Synthesize modulator
    if fm == 0.0
        modulator = zeros(length(carrier))
    else
        modulator = pure_tone(fm, ϕ, dur, fs)
    end

    # Combine carrier and modulator and scale
    stim = (1.0 .+ m .* modulator) .* carrier
    scale_dbspl(stim, level_pre_mod)
end