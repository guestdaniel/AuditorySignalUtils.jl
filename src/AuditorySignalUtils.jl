module AuditorySignalUtils

using Statistics

export amplify, amplify!, cosine_ramp, cosine_ramp!, dbspl, LogRange, pure_tone, rms, 
       scale_dbspl, scale_dbspl!

"""
    amplify(signal, dB)

Amplifies a signal (in terms of power) by an amount in dB (or attenuates the signal if dB\
 is negative)
"""
function amplify(signal::AbstractVector{T}, dB::T)::AbstractVector{T} where {T<:Real}
    signal .* 10^(dB/20)
end


"""
    amplify!(signal, dB)

Amplifies a signal (in terms of power) by an amount in dB (or attenuates the signal if dB is negative).

Note, this version of the function operates in-place (i.e., no extra memory is allocated and the input is modified in-place)
"""
function amplify!(signal::AbstractVector{T}, dB::T)::AbstractVector{T} where {T<:Real}
    signal .*= 10^(dB/20)
    return signal
end


"""
    cosine_ramp(signal, dur_ramp, fs)

Applies a raised-cosine ramp to the input signal, where dur_ramp is the duration of each
samp segment in seconds
"""
function cosine_ramp(signal::AbstractVector{T}, dur_ramp::T, fs::T)::AbstractVector{T} where {T<:Real}
    t = LinRange(0, 0.25, Int(floor(fs*dur_ramp)))  # f = 1 Hz, dur = 1/4 cycle, 0 -> 1
    ramp_segment = sin.(2*pi*t).^2
    ramp = [ramp_segment; ones(length(signal) - length(ramp_segment)*2); reverse(ramp_segment)]
    return signal .* ramp
end


"""
    cosine_ramp!(signal, dur_ramp, fs)

Applies a raised-cosine ramp to the input signal, where dur_ramp is the duration of each
samp segment in seconds
"""
function cosine_ramp!(signal::AbstractVector{T}, dur_ramp::T, fs::T)::AbstractVector{T} where {T<:Real}
    t = LinRange(0, 0.25, Int(floor(fs*dur_ramp)))  # f = 1 Hz, dur = 1/4 cycle, 0 -> 1
    ramp_segment = sin.(2*pi*t).^2
    ramp = [ramp_segment; ones(length(signal) - length(ramp_segment)*2); reverse(ramp_segment)]
    signal .*= ramp
    return signal
end


"""
    dbspl(signal)::AbstractFloat

Calculates the dB SPL value of a signal (assuming that the units of the signal are in Pa)
"""
function dbspl(signal::AbstractVector{T})::T where {T<:Real}
    20*log10(rms(signal)/20e-6)
end


"""
    LogRange(a, b, n)::Range
"""
function LogRange(a::T, b::T, n::Int) where {T<:Real}
    10 .^ LinRange(log10(a), log10(b), n)
end


"""
    pure_tone(freq, phase, dur, fs)

Synthesizes a pure tone with specified frequency, phase offset (in radians), duration, and sampling rate.
"""
function pure_tone(freq::T, phase::T, dur::T, fs::T)::AbstractVector{T} where {T<:Real}
    sin.(2*pi .* freq .* (0.0:(1/fs):prevfloat(dur)) .+ phase)
end


"""
    rms(signal)::AbstractFloat

Calculates the room-mean-square (RMS) of a signal
"""
function rms(signal::AbstractVector{T})::T where {T<:Real}
    sqrt(mean(signal.^2))
end


"""
    scale_dbspl(signal, dB)

Adjusts a signals level to be a certain level in dB SPL
"""
function scale_dbspl(signal::AbstractVector{T}, dB::T)::AbstractVector{T} where {T<:Real}
    curr_dB = dbspl(signal)
    delta_dB = dB - curr_dB
    amplify(signal, delta_dB)
end

"""
    scale_dbspl!(signal, dB)

Adjusts a signals level to be a certain level in dB SPL
"""
function scale_dbspl!(signal::AbstractVector{T}, dB::T)::AbstractVector{T} where {T<:Real}
    amplify!(signal, dB - dbspl(signal))
end


end # module
