export klatt1980, klatt1980_parallel

"""
    klatt_resonator_coefs(ff, bw, fs)

Compute coefficients for Klatt (1980) resonator given ff, bw, and fs
"""
function klatt_resonator_coefs(ff, bw, fs)
    C = -exp(-2π * bw/fs)
    B = 2 * exp(-π * bw/fs) * cos(2 * π * ff/fs)
    A = 1 - B - C
    b = [A]
    a = [1.0, -B, -C]
    return b, a
end

"""
    klatt_antiresonator_coefs(ff, bw, fs)

Compute coefficients for Klatt (1980) antiresonator given ff, bw, and fs
"""
function klatt_antiresonator_coefs(ff, bw, fs)
    C = -exp(-2π * bw/fs)
    B = 2 * exp(-π * bw/fs) * cos(2 * π * ff/fs)
    A = 1 - B - C
    A′ = 1.0/A
    B′ = -B/A
    C′ = -C/A
    b = [A′, -B', -C']
    a = [1.0]
    return b, a
end

"""
    klatt_resonator(x, ff, bw, fs)

Filter signal x with Klatt (1980) resonator given parameters (see `klatt_resonator_coefs` for details)
"""
klatt_resonator(x, args...) = klatt_resonator!(zeros(size(x)), x, args...)
klatt_resonator!(y, x, args...) = filt!(y, klatt_resonator_coefs(args...)..., x)

"""
    klatt_antiresonator(x, ff, bw, fs)

Filter signal x with Klatt (1980) antiresonator given parameters (see `klatt_antiresonator_coefs` for details)
"""
klatt_antiresonator(x, args...) = klatt_antiresonator!(zeros(size(x)), x, args...)
klatt_antiresonator!(y, x, args...) = filt!(y, klatt_antiresonator_coefs(args...)..., x)

"""
    klatt1980(; kwargs...)

Synthesizes a vowel according to a simplified cascade architecture from Klatt (1980).

Synthesizes a Klatt vowel according to a modified version of Klatt (1980). First, a
perfectly periodic pulse train with F0 is synthesized. This pulse train is optionally
filtered by a glottal pole and zero (RGP and RGZ, respectively, in Klatt [1980]), 
an arbitrary number of formants in cascade, and then the radiation characteristic associated
with the lips.
"""
function klatt1980(; 
    dur=1.0, 
    f0=100.0, 
    ff=[300.0, 1000.0], 
    bw=fill(100.0, size(ff)), 
    fs=100e3, 
    rgp=true,
    rgz=true,
    rad=true,
)
    # Construct stimulus impulse train
    stim = zeros(samples(dur, fs))
    stim[1:samples(1/f0, fs):end] .= 1.0

    # Apply RGP + RGZ
    if rgp stim = klatt_resonator(stim, 0.0, 100.0, fs) end
    if rgz stim = klatt_antiresonator(stim, 1500.0, 6000.0, fs) end

    # Apply each filter in the cascade path
    for (_ff, _bw) in zip(ff, bw)
        stim = klatt_resonator(stim, _ff, _bw, fs)
    end

    # Filter stimulus again to reflect voicing characteristic
    if rad stim = filt([1.0, -1.0], [1.0], stim) end

   return stim
end

"""
    klatt1980_parallel(; kwargs...)

Synthesizes a vowel according to a simplified parallel architecture from Klatt (1980).

Synthesizes a Klatt vowel according to a modified version of Klatt (1980). First, a
perfectly periodic pulse train with F0 is synthesized. This pulse train is optionally
filtered by a glottal pole and zero (RGP and RGZ, respectively, in Klatt [1980]), 
an arbitrary number of formants *in parallel*, and then the radiation characteristic associated
with the lips.
"""
function klatt1980_parallel(; 
    dur=1.0, 
    f0=100.0, 
    ff=[300.0, 1000.0], 
    bw=fill(100.0, size(ff)), 
    autoflatten=false,
    amp=autoflatten ? 10 .^ (-6.0 .* (log2.(ff ./ ff[1])) ./ 20.0) : ones(size(ff)),
    fs=100e3, 
    rgp=true,
    rgz=true,
    rad=true,
)
    # Construct stimulus impulse train
    stim = zeros(samples(dur, fs))
    stim[1:samples(1/f0, fs):end] .= 1.0

    # Apply RGP + RGZ
    if rgp stim = klatt_resonator(stim, 0.0, 100.0, fs) end
    if rgz stim = klatt_antiresonator(stim, 1500.0, 6000.0, fs) end

    # Synthesize each waveform in the parallel path, then combine
    if length(ff) == 0
        output = stim
    else
        output = zeros(samples(dur, fs))
        for (_ff, _bw, _amp) in zip(ff, bw, amp)
            output .+= _amp .* klatt_resonator(stim, _ff, _bw, fs)
        end
    end

    # Filter stimulus again to reflect voicing characteristic
    if rad output = filt([1.0, -1.0], [1.0], output) end
    return output
end
