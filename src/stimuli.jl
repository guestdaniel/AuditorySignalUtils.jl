export pure_tone, gaussian_noise, gaussian_noise_sl, sam_noise, sam_noise_sl

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