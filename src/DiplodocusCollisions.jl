module DiplodocusCollisions

export BinaryInteractionIntegration, BinaryFileLoad_All, DoesConserve, BinaryFileLoad_Matrix
export EmissionInteractionIntegration, EmissionFileLoad_Matrix, EmissionFileLoad_All
export UserBinaryParameters, UserEmissionParameters

    using JLD2
    using Base.Threads
    using BenchmarkTools
    using Bessels
    using ProgressMeter

    # include Common files
        include("Common/Constants.jl")
        include("Common/RandomPoints.jl")
        include("Common/Location.jl")
        include("Common/UserParameters.jl")

    # include Binary files
        #include("Binary/Structs.jl")
        include("Binary/Arrays.jl")
        include("Binary/Averaging.jl")
        include("Binary/DifferentialCrossSectionFunctions.jl")
        include("Binary/MomentumValues.jl")
        include("Binary/MandelstramChecks.jl")
        include("Binary/GainLossValue.jl")
        include("Binary/UsefulGridValueFunctions.jl")
        include("Binary/MomentumSpaceFactors.jl")
        include("Binary/Sampling.jl")
        include("Binary/BinaryInteractionIntegration.jl")
        include("Binary/BinaryMonteCarlo.jl")       
        include("Binary/BinaryMonteCarlo_Debug.jl")   
        include("Binary/BinaryDataReading.jl")
        
    # include Synchrotron functions
        include("Emission/Arrays.jl")
        include("Emission/Averaging.jl")
        include("Emission/EmissionKernels.jl")
        include("Emission/MomentumSpaceFactors.jl")
        include("Emission/Sampling.jl")
        include("Emission/EmissionInteractionIntegration.jl")
        include("Emission/EmissionMonteCarlo.jl")
        include("Emission/EmissionMonteCarlo_Debug.jl")
        include("Emission/EmissionDataReading.jl")

end

