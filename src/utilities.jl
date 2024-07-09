# Amplification, level setting, level measurements, etc.
export amplify, amplify!, dbspl, scale_dbspl, scale_dbspl!, sl_to_ol

# Signal ramping
export cosine_ramp, cosine_ramp!

# Log ranges, octave ranges, calculating octaves, etc.
export LogRange, OctRange, octs, logtimerange

# Zero-padding and other generic utilities
export zero_pad, zero_pad!, silence, withisi

# Conversion between samples, times, durations, etc.
export samples, sampleat, timeat, timevec

# Compute octaves
octs(freq, shift) = freq * 2.0 ^ shift

# Convert duration to samples 
samples(time, fs) = Int(round(time*fs))

# Convert sample index to time  index
timeat(sample, fs) = (sample-1)/fs

# Convert time index to sample index
sampleat(time, fs) = Int(round(time*fs + 1))

# Create time vector
timevec(dur::Float64, fs) = (0.0:(1/fs):nextfloat(dur-1/fs))
timevec(samples::Int64, fs) = (0.0:(1/fs):nextfloat(samples/fs-1/fs))
timevec(x::Vector, fs) = (0.0:(1/fs):nextfloat(length(x)/fs-1/fs))

# Create range of times on log2 scale
logtimerange(time, low=-3, high=3, spacing=1) = time .* 2.0 .^ (low:spacing:high)

# Create period of silence
silence(time, fs) = zeros(samples(time, fs))

# Put silence between two vectors
withisi(x, y; isi=0.10, fs=100e3) = vcat(x, silence(isi, fs), y)

"""
    amplify(x, dB)

Amplifies (or attenuates) signal by factor of 10^(`dB`/20)
"""
amplify!(x, dB) = x .*= 10.0^(dB/20.0)
amplify(x, dB) = x .* 10.0^(dB/20.0)

"""
    cosine_ramp(x, dur_ramp, fs)

Applies raised-cosine ramp of dur `dur_ramp` (s) to input signal of sampling rate `fs`
"""
function cosine_ramp(x, dur_ramp, fs)
    t = LinRange(0, prevfloat(0.25), Int(floor(fs*dur_ramp)))  # f = 1 Hz, dur = 1/4 cycle, 0 -> 1
    ramp_segment = sin.(2π .* t).^2
    ramp = [ramp_segment; ones(length(x) - length(ramp_segment)*2); reverse(ramp_segment)]
    return x .* ramp
end

function cosine_ramp!(x, dur_ramp, fs)
    t = LinRange(0, prevfloat(0.25), Int(floor(fs*dur_ramp)))  # f = 1 Hz, dur = 1/4 cycle, 0 -> 1
    ramp_segment = sin.(2π .* t).^2
    ramp = [ramp_segment; ones(length(x) - length(ramp_segment)*2); reverse(ramp_segment)]
    x .*= ramp
    return x
end

"""
    dbspl(x)

Calculate signal's RMS level in dB SPL re: 20 μPa
"""
dbspl(x) = 20.0*log10(DSP.rms(x)/20e-6)

"""
    LogRange(a, b, n)

Creates a vector with n elements spaced logarithmically from a to b
"""
function LogRange(a::T, b::T, n::Int) where {T<:Real}
    if n == 1
        exp(1/2 * (log(a) + log(b)))
    else
        exp.(LinRange(log(a), log(b), n))
    end
end

"""
    OctRange(a, b, n)

Creates a vector with n elements spaced from `l` octs below to `u` octs above `x`
"""
function OctRange(x::T, l::T, u::T, n::Int) where {T<:Real}
    if n == 1
        error("Can't octrange one value")
    else
        LogRange(octs(x, l), octs(x, u), n)
    end
end

"""
    scale_dbspl(x, level)

Adjusts signal's level to be `level` in dB SPL re: 20μPa
"""
scale_dbspl(x, level) = amplify(x, isnan(level - dbspl(x)) ? -Inf : level - dbspl(x))
scale_dbspl!(x, level) = amplify!(x, isnan(level - dbspl(x)) ? -Inf : level - dbspl(x))

"""
    sl_to_ol(level, fs)

Converts broadband spectrum `level` with sampling rate `fs` to overall level
"""
sl_to_ol(level, fs) = level + 10*log10(fs/2)