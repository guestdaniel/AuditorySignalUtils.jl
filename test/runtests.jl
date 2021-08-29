using Test
using AuditorySignalUtils

# Declare various constants that hold across all tests
fs = 5000.0
dur = 1.0
precision = -3

# Define a function to check approximate equality
function approx_equal(x1, x2, error_size)
    abs(x1 - x2) <= 10.0^error_size
end

# Define a shorter function that handles approximate equality with a constant precision across tests
function ae(x1, x2)
    approx_equal(x1, x2, precision)
end

# Test the approx_equal function
@test approx_equal(0, 0.0000000001, -5)
@test approx_equal(0, 1, -5) == false

# Test that pure_tone() properly synthesizes tones with phase shifts in radians (0, pi/2, pi, 3pi/2) work correctly
@testset "pure tones" begin
    @test ae(AuditorySignalUtils.pure_tone(5.0, 0.0, dur, fs)[1], 0)
    @test ae(AuditorySignalUtils.pure_tone(5.0, pi/2, dur, fs)[1], 1)
    @test ae(AuditorySignalUtils.pure_tone(5.0, pi, dur, fs)[1], 0)
    @test ae(AuditorySignalUtils.pure_tone(5.0, 3*pi/2, dur, fs)[1], -1)
end

# Test that rms() correctly returns the rms values of an unscaled sinusoid
@test ae(AuditorySignalUtils.rms(AuditorySignalUtils.pure_tone(5.0, 0.0, dur, fs)), 0.707)

# Test that dbspl() correctly correctly calculates levels of reference stimuli
@test ae(AuditorySignalUtils.dbspl(AuditorySignalUtils.pure_tone(5.0, 0.0, dur, fs)), 90.968)

# Test that amplify() produces correct changes in power and amplitude of sinusoid
@testset "amplification" begin
    # Test amplification
    @test begin
        # Synthesize pure tone and then increment it by 10 dB
	      baseline_signal = AuditorySignalUtils.pure_tone(5.0, 0.0, dur, fs)
        amplified_signal = AuditorySignalUtils.amplify(baseline_signal, 6.0)
        # Check that amplitude and power ratios are correct
        ae(maximum(amplified_signal)/maximum(baseline_signal), 1.995) && ae(maximum(amplified_signal)^2/maximum(baseline_signal)^2, 3.981)
    end
    # Test attenuation
    @test begin
        # Synthesize pure tone and then increment it by 10 dB
	      baseline_signal = AuditorySignalUtils.pure_tone(5.0, 0.0, dur, fs)
        amplified_signal = AuditorySignalUtils.amplify(baseline_signal, -6.0)
        # Check that amplitude and power ratios are correct
        ae(maximum(amplified_signal)/maximum(baseline_signal), 0.501) && ae(maximum(amplified_signal)^2/maximum(baseline_signal)^2, 0.251)
    end
end
