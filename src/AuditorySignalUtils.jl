module AuditorySignalUtils
export amplify, cosine_ramp, dbspl, LogRange, pure_tone, rms, scale_dbspl 

using Statistics

"""
    amplify(signal, dB)

Amplifies a signal's power by a certain amount in dB (or attenuates the signal if dB is negative)
"""
function amplify(signal::Array{<:AbstractFloat, 1}, dB)::Array{<:AbstractFloat, 1}
    signal .* 10^(dB/20)
end


"""
    cosine_ramp(signal, dur_ramp, fs)

Applies a raised-cosine ramp to the input signal, where dur_ramp is the duration of each
samp segment in seconds
"""
function cosine_ramp(signal::Array{<:AbstractFloat, 1}, dur_ramp, fs)::Array{<:AbstractFloat, 1}
    t = LinRange(0, 0.25, Int(floor(fs*dur_ramp)))  # f = 1 Hz, dur = 1/4 cycle, 0 -> 1
    ramp_segment = sin.(2*pi*t).^2
    ramp = [ramp_segment; ones(length(signal) - length(ramp_segment)*2); reverse(ramp_segment)]
    return signal .* ramp
end


"""
    dbspl(signal)::AbstractFloat

Calculates the dB SPL value of a signal (assuming that the units of the signal are in Pascals)
"""
function dbspl(signal::Array{<:AbstractFloat, 1})::AbstractFloat
    20*log10(rms(signal)/20e-6)
end


"""
    LogRange(a, b, n)::Range
"""
function LogRange(a, b, n::Int)
    10 .^ LinRange(log10(a), log10(b), n)
end


"""
    pure_tone(freq, phase, dur, fs)

Synthesizes a pure tone with specified frequency, phase offset (in radians), duration, and sampling rate.
"""
function pure_tone(freq, phase, dur, fs)::Array{<:AbstractFloat, 1}
    sin.(2*pi*freq*LinRange(0, dur, Int64(dur*fs)) .+ phase)
end


"""
    rms(signal)::AbstractFloat

Calculates the room-mean-square (RMS) of a signal
"""
function rms(signal::Array{<:AbstractFloat, 1})::AbstractFloat
    sqrt(mean(signal.^2))
end


"""
    scale_dbspl(signal, dB)

Adjusts a signals level to be a certain level in dB SPL
"""
function scale_dbspl(signal::Array{<:AbstractFloat, 1}, dB)::Array{<:AbstractFloat, 1}
    curr_dB = dbspl(signal)
    delta_dB = dB - curr_dB
    amplify(signal, delta_dB)
end

end # module
