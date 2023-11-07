module AuditorySignalUtils

using Statistics
using DSP: rms

export amplify, amplify!, 
       cosine_ramp, cosine_ramp!, 
       dbspl, 
       LogRange, 
       pure_tone, 
       scale_dbspl, scale_dbspl!, 
       zero_pad, zero_pad!

samples(time, fs) = Int(round(time, fs))
timeat(sample, fs) = sample*fs
time(dur, fs) = (0.0:(1/fs):(dur-1/fs))

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

Calculate signal's RMS level in dB SPL re: 20μPa
"""
dbspl(x) = 20.0*log10(rms(x)/20e-6)

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
    pure_tone(freq, phase, dur, fs)

Synthesize pure tone; `freq` (Hz), phase `ϕ` (rads), `dur` (s), sample rate `fs` (Hz)
"""
pure_tone(freq, ϕ, dur, fs) = sin.(2π .* freq .* time(dur, fs) .+ ϕ)

"""
    scale_dbspl(x, level)

Adjusts signal's level to be `level` in dB SPL re: 20μPa
"""
scale_dbspl(x, level) = amplify(x, level - dbspl(x))
scale_dbspl!(x, level) = amplify!(x, level - dbspl(x))

end # module


