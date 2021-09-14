module AuditorySignalUtils
using Statistics

"""
    amplify(signal, dB)

Amplifies a signal's power by a certain amount in dB (or attenuates the signal if dB is negative)
"""
function amplify(signal::Array{Float64, 1}, dB::Float64)::Array{Float64, 1}
    signal .* 10^(dB/20)
end


"""
    dbspl(signal)::Float64

Calculates the dB SPL value of a signal (assuming that the units of the signal are in Pascals)
"""
function dbspl(signal::Array{Float64, 1})::Float64
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
function pure_tone(freq::Float64, phase::Float64, dur::Float64, fs::Float64)::Array{Float64, 1}
    sin.(2*pi*freq*LinRange(0, dur, Int64(dur*fs)) .+ phase)
end

function pure_tone(freq::Float64, phase::Irrational, dur::Float64, fs::Float64)::Array{Float64, 1}
    sin.(2*pi*freq*LinRange(0, dur, Int64(dur*fs)) .+ phase)
end


"""
    rms(signal)::Float64

Calculates the room-mean-square (RMS) of a signal
"""
function rms(signal::Array{Float64, 1})::Float64
    sqrt(mean(signal.^2))
end


"""
    scale_dbspl(signal, dB)

Adjusts a signals level to be a certain level in dB SPL
"""
function scale_dbspl(signal::Array{Float64, 1}, dB::Float64)::Array{Float64, 1}
    curr_dB = dbspl(signal)
    delta_dB = dB - curr_dB
    amplify(signal, delta_dB)
end

end # module
