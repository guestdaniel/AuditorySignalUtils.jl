module AuditorySignalUtils

using Statistics
using DSP: rms

export amplify, amplify!, cosine_ramp, cosine_ramp!, dbspl, LogRange, pure_tone, 
       scale_dbspl, scale_dbspl!, zero_pad

"""
    amplify(x, dB)

Amplifies (or attenuates) a signal by an amount in dB

Changes the power in signal `x` by decibel factor `dB`. When `dB` is positive, the power in
the signal is increased (i.e., signal is amplified) and when `dB` is negative, the power in
the signal is decreased (i.e, signal is attenuated). Output signal and input signal have an
intensity ratio of 10^(`dB`/10) and an amplitude ratio of 10^(`dB`/20)
"""
function amplify(x::AbstractVector{T}, dB::T)::AbstractVector{T} where {T<:AbstractFloat}
    x .* 10.0^(dB/20.0)
end

function amplify!(x::AbstractVector{T}, dB::T)::AbstractVector{T} where {T<:AbstractFloat}
    x .*= 10.0^(dB/20.0)
    return x
end

"""
    cosine_ramp(x, dur_ramp, fs)

Applies a raised-cosine ramp to an input signal

Applies a raised-cosine ramp to the input signal `x` where the ramp has duration `dur_ramp`
(s) and sampling rate `fs` (Hz)
"""
function cosine_ramp(x::AbstractVector{T}, dur_ramp::T, fs::T)::AbstractVector{T} where {T<:Real}
    t = LinRange(0, prevfloat(0.25), Int(floor(fs*dur_ramp)))  # f = 1 Hz, dur = 1/4 cycle, 0 -> 1
    ramp_segment = sin.(2π .* t).^2
    ramp = [ramp_segment; ones(length(x) - length(ramp_segment)*2); reverse(ramp_segment)]
    return x .* ramp
end

function cosine_ramp!(x::AbstractVector{T}, dur_ramp::T, fs::T)::AbstractVector{T} where {T<:Real}
    t = LinRange(0, prevfloat(0.25), Int(floor(fs*dur_ramp)))  # f = 1 Hz, dur = 1/4 cycle, 0 -> 1
    ramp_segment = sin.(2π .* t).^2
    ramp = [ramp_segment; ones(length(x) - length(ramp_segment)*2); reverse(ramp_segment)]
    x .*= ramp
    return x
end

"""
    dbspl(x)

Calculates the dB SPL value of a signal

Calculates the dB SPL of the whole signal `x` assuming that the signals unit's are Pa
"""
function dbspl(x::AbstractVector{T})::T where {T<:Real}
    20*log10(rms(x)/20e-6)
end

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

Synthesize a pure tone 

Synthesize pure tone with with specified frequency `f` (Hz), phase offset `ϕ` (rads),
duration `dur` (s), and sampling rate `fs` (Hz) 
"""
function pure_tone(freq::T, ϕ::T, dur::T, fs::T)::AbstractVector{T} where {T<:Real}
    sin.(2π .* freq .* (0.0:(1/fs):(dur-1/fs)) .+ ϕ)
end

"""
    pure_tone(freq, phase, dur, fs)

Synthesize a pure tone 

Synthesize pure tone with with specified frequency `f` (Hz), phase offset `ϕ` (rads),
duration `dur` (s), and sampling rate `fs` (Hz) 
"""
function pure_tone(freq::T, ϕ::M, dur::T, fs::T)::AbstractVector{T} where {T<:Real, M<:Irrational}
    sin.(2π .* freq .* (0.0:(1/fs):(dur-1/fs)) .+ ϕ)
end

"""
    scale_dbspl(x, dB)

Adjusts a signal's level to be a certain level in dB SPL

Modified signal `x` to have level `dB` in dB SPL 
"""
function scale_dbspl(x::AbstractVector{T}, dB::T)::AbstractVector{T} where {T<:Real}
    curr_dB = dbspl(x)
    delta_dB = dB - curr_dB
    amplify(x, delta_dB)
end

function scale_dbspl!(x::AbstractVector{T}, dB::T)::AbstractVector{T} where {T<:Real}
    amplify!(x, dB - dbspl(x))
end

end # module

