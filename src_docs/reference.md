# Reference

Here, you can see the documentation for each function provided by AuditorySignalUtils.jl.

```@meta
CurrentModule = AuditorySignalUtils
```

## Synthesis

Synthesize stimuli. 
Right now just ramps and pure tones.

```@docs
pure_tone
cosine_ramp
```

## Digital signal processing tools

Various functions to handle things like calculating rms, setting levels, etc.

```@docs
rms
dbspl
scale_dbspl
amplify
```

A few functions include in-place versions to optimize performance. 
These versions are marked with a trailing exclamation mark (!) per standard Julia style.

```@docs
scale_dbspl!
amplify!
```

## Other utilities

Various functions I frequently need and don't want to redefine in every file.

```@docs
LogRange
```
