# Amplification, level setting, level measurements, etc.
export amplify, amplify!, dbspl, scale_dbspl, scale_dbspl!, sl_to_ol

# Signal ramping
export cosine_ramp, cosine_ramp!

# Log ranges, octave ranges, calculating octaves, etc.
export LogRange, OctRange, octs, logtimerange

# Zero-padding and other generic utilities
export zero_pad, zero_pad!, silence, withisi, mixat

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

# Convert start and stop time to corresponding sample range, under the normal convention
# that stop time is exclusive w.r.t to the nearest sample, i.e., stop time is really 
# stop time - 1/fs
function samples(start, stop, fs)
    start = sampleat(start, fs)
    stop = samples(stop, fs)
    if start > stop
        error("Start time must be less than stop time")
    end
    return start:stop
end

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

# Mix vector y with vector x starting at d seconds after the start of x
function mixat(x, y, d; fs=100e3)
    idx = sampleat(d, fs)
    z = deepcopy(x)
    z[idx:(idx+length(y)-1)] .+= y
    z
end

# Basic signal processing utilitiesj
hann(n, N) = (1/2) * (1 - cos(2π*n/N))

"""
    amplify(x, dB)

Amplifies (or attenuates) signal by factor of 10^(`dB`/20)
"""
amplify!(x, dB) = x .*= 10.0^(dB/20.0)
amplify(x, dB) = x .* 10.0^(dB/20.0)

"""
    cosine_ramp(x, dur_ramp, fs)

Applies raised-cosine ramp of dur `dur_ramp` (s) to input signal of sampling rate `fs` (Hz)

Applies a raised-cosine ramp to the input `x` with specified duration and sampling rate.
Note that this ramp is the square of a cosine ramp, or equivalently the Hann function:
    1/2 * [1 - cos(2πn/N)], 0 ≤ n ≤ N
evaluated from 0 to N where N is twice the length of the ramp in samples:
"""
function cosine_ramp(x, dur_ramp, fs)
    y = deepcopy(x)
    cosine_ramp!(y, dur_ramp, fs)
end

function cosine_ramp!(x, dur_ramp, fs)
    len = samples(dur_ramp, fs)
    r = hann.(0:len, len*2)
    r = vcat(r, ones(length(x) - 2*length(r)), reverse(r))
    x .*= r
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