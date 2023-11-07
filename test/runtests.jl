using Test
using AuditorySignalUtils
using DSP: rms

# Declare various constants that hold across all tests
const fs = 100e3

# Test: amplify correctly adjusts power and amplitude ratios
@testset "amplify" begin
    for dB in -50.0:5.0:50.0
        x_ref = pure_tone(1000.0, 0.0, 0.1, fs)
        x_amp = amplify(x_ref, dB)
        amplitude_ratio = rms(x_amp) / rms(x_ref)
        power_ratio = rms(x_amp)^2 / rms(x_ref)^2
        @test (
            (amplitude_ratio ≈ 10^(dB/20)) & 
            (power_ratio ≈ 10^(dB/10))
        )
    end
end

# Test: `pure_tone` properly phase shifts a sine tone as function of input phase 
@testset "Pure tone synthesis" begin
    for ϕ in LinRange(0.0, 2π, 20)
        @test pure_tone(5.0, ϕ, 0.01, fs)[1] ≈ sin(ϕ)
    end
end

# Test: `dbspl` returns the correct value for various reference sounds
@testset "dB SPL" begin
    @test dbspl(zeros(5000)) == -Inf
    @test dbspl(pure_tone(1000.0, 0.0, 1.0, fs)) ≈ 90.969 atol=0.01
    @test begin
        x = pure_tone(1000.0, 0.0, 1.0, fs)
        x = x ./ rms(x)
        ≈(dbspl(x), 93.979; atol=0.01)
    end
end

