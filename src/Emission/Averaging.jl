"""
    WeightedAverageGainEmission!(GainMatrix2,OldGainMatrix2,GainTally2,OldGainTally2,GainMatrix3,OldGainMatrix3,GainTally3,OldGainTally3)

Computes the integral estimate by weighted average of the old and new gain matrices. Mutating the old gain and tally terms.
"""
function WeightedAverageGainEmission!(GainMatrix2::Array{Float64,6},OldGainMatrix2::Array{Float64,6},GainTally2::Array{UInt32,6},OldGainTally2::Array{UInt32,6},GainMatrix3::Array{Float64,6},OldGainMatrix3::Array{Float64,6},GainTally3::Array{UInt32,6},OldGainTally3::Array{UInt32,6})

    # weighted average
    @. OldGainMatrix2 = (GainMatrix2*GainTally2+OldGainMatrix2*OldGainTally2)/(GainTally2+OldGainTally2) 
    @. OldGainMatrix3 = (GainMatrix3*GainTally3+OldGainMatrix3*OldGainTally3)/(GainTally3+OldGainTally3)

    replace!(OldGainMatrix2,NaN=>0e0)
    replace!(OldGainMatrix3,NaN=>0e0)

    # adding tallies 
    @. OldGainTally2 += GainTally2
    @. OldGainTally3 += GainTally3

end

"""
    WeightedAverageLossEmission!(LossMatrix1,OldLossMatrix1,LossTally1,OldLossTally1)

Computes the integral estimate by weighted average of the old and new gain matrices. Mutating the old gain and tally terms.
"""
function WeightedAverageLossEmission!(LossMatrix1::Array{Float64,3},OldLossMatrix1::Array{Float64,3},LossTally1::Array{UInt32,3},OldLossTally1::Array{UInt32,3})

    # weighted average 
    @. OldLossMatrix1 = (LossMatrix1*LossTally1+OldLossMatrix1*OldLossTally1)/(LossTally1+OldLossTally1)

    replace!(OldLossMatrix1,NaN=>0e0)

    # adding tallies
    @. OldLossTally1 += LossTally1

end