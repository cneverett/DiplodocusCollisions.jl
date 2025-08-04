<img align="left" width="75" height="75" src="https://github.com/cneverett/BoltzmannCollisionIntegral.jl/blob/main/docs/src/assets/logo.svg" alt="BoltzmannCollisionIntegral.jl icon">

# BoltzmannCollisionIntegral.jl

[![][docs-latest-img]][docs-latest-url]
[![Build Status][gha-img]][gha-url]

[docs-latest-img]: https://img.shields.io/badge/Docs-Stable-lightgrey.svg
[docs-latest-url]: https://cneverett.github.io/BoltzmannCollisionIntegral.jl/dev/

[gha-img]: https://github.com/cneverett/BoltzmannCollisionIntegral.jl/actions/workflows/CI.yml/badge.svg?branch=main
[gha-url]: https://github.com/cneverett/BoltzmannCollisionIntegral.jl/actions/workflows/CI.yml?query=branch%3Amain

> [!CAUTION]
> This project is undergoing an major revision to align with the Diplodocus package framework



## Overview

`BoltzmannCollisionIntegral.jl` is a [Julia](http://julialang.org/) package for the evaluation of the relativistic Boltzmann collision integral for binary interactions $(12\rightleftharpoons34)$ (see [Everett+Cotter 2024](https://doi.org/10.1093/rasti/rzae036)):
```math
    C(\boldsymbol{p}_3)=\int\frac{\mathrm{d}^3\boldsymbol{p}_1}{p_1^0}\frac{\mathrm{d}^3\boldsymbol{p}_2}{p_2^0}\frac{\mathrm{d}^3\boldsymbol{p}_4}{p_4^0}\left[\frac{f(\boldsymbol{p}_1)f(\boldsymbol{p}_2)}{1+\delta_{12}}W(p_1^\mu,p_2^\mu|p_3^\mu,p_4^\mu)- \frac{f(\boldsymbol{p}_3)f(\boldsymbol{p}_4)}{1+\delta_{12}}W(p_3^\mu,p_4^\mu|p_1^\mu,p_2^\mu)\right],
```
via momentum discretisation and Monte-Carlo sampling. The distribution functions $f(\boldsymbol{p})$ for the particles involved are assumed to be anisotropic (only axisymmetric is currently implemented). The collision integral is split into two components, an emission spectra dictating the rate of gain of particles emerging from a given binary interaction and an absorption spectra dictating the rate of loss of particles from the same interaction.

Evaluation is performed by assuming that the distribution function is constant over some set of discrete domains in phase space and then integrating over those domains using a Monte-Carlo method. For more information see the [Documentation](https://cneverett.github.io/BoltzmannCollisionIntegral.jl/).

Evaluation currently supports both single and multithreaded operation (with multi-CPU acceleration planned). Data is exported in the [JLD2](https://github.com/JuliaIO/JLD2.jl) file format.

Since v2.0.0 the evaluation of synchrotron emissions from charged particles orbiting a uniform magnetic field directed along the symmetry axis is also supported. This is achieved via the integration of the following:
```math
\left(\frac{\partial f(\boldsymbol{p}_1)}{\partial t}\right)_{Sync}=\frac{Z}{B}\frac{3c^5\sigma_Tm_e^5}{4\pi\hbar^3\mu_0e^3}\int\mathrm{d}^3\boldsymbol{p}_2\frac{1}{p_1p^0_2}\left[\left(\frac{p^0_2\cos\theta_1-p_2\cos\theta_2}{\sin\theta_1}\right)^2J_{\omega/\omega_0}(x)^2+(p_2\sin\theta_2)^2 J'_{\omega/\omega_0}(x)^2\right]f(\boldsymbol{p}_2).
```

## Usage
`BoltzmannCollisionIntegral.jl` is available to download from the [Julia package
manager](https://pkgdocs.julialang.org/v1/). Inside a Julia session, enter the package manager with `]`, then run the command 
```julia 
pkg> add BoltzmannCollisionIntegral
```
finally load the package by running
```julia
using BoltzmannCollisionIntegral
```

To perform an evaluation of the emission and absorption spectra, an example script `Run_Integration.jl` for selecting the binary interaction, discrete phase space bounds and integration conditions is located under the `src/Common/` folder of the package. It is recommended to copy this script and place it in your working folder and edit the fields as required. Then simply run
```julia-repl
include("Run_Integration.jl")
``` 

To perform an evaluation of the synchrotron emission spectra, an example script `Run_Integration_Sync.jl` for selecting the emitting particle, discrete phase space bounds and integration conditions is located under the `src/Synchrotron/Common/` folder of the package. It is recommended to copy this script and place it in your working folder and edit the fields as required. Then simply run
```julia-repl
include("Run_Integration_Sync.jl")
``` 

See [Getting Started](https://cneverett.github.io/BoltzmannCollisionIntegral.jl/dev/quickstart/) for in depth detail.