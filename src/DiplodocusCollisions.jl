module DiplodocusCollisions

export BinaryInteractionIntegration, BinaryFileLoad_All, DoesConserve, BinaryFileLoad_Matrix #, BinaryFileLoad_Matrix_ISO
export SyncEvaluateSerial, SpectraEvaluateMultiThreadEmission, fload_All_Sync, fload_Matrix_Sync, fload_Matrix_SyncISO
export UserBinaryParameters 

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
        include("Binary/DataReading.jl")
        
    # include Synchrotron functions
        include("Emission/Common/Arrays.jl")
        include("Emission/Common/Averaging.jl")
        include("Emission/Common/SynchrotronKernel.jl")
        include("Emission/Common/MomentumSpaceFactors.jl")
        include("Emission/Common/Sampling.jl")
        # include serial methods
        include("Emission/Serial/SyncMonteCarlo_Serial.jl")
        include("Emission/Serial/SyncIntegration_Serial.jl")
        # include parallel methods
        include("Emission/MultiThread/SyncMonteCarlo_MultiThread.jl")
        include("Emission/MultiThread/SyncIntegration_MultiThread.jl")
        # include data reading functions for export
        include("Emission/Common/SyncDataReading.jl")

end

