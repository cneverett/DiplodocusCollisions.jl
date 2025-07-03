"""
    WeightedAverageGainBinary!(GainMatrix3,OldGainMatrix3,GainTally3_K,OldGainTally3_K,GainMatrix4,OldGainMatrix4,GainTally4_K,OldGainTally4_K)

Computes the integral estimate by weighted average of the old and new gain matrices. Mutating the old gain and tally terms.
"""
function WeightedAverageGainBinary!(GainMatrix3::Array{Float64,9},OldGainMatrix3::Array{Float64,9},GainTally3_K::AbstractArray{UInt32,9},OldGainTally3::Array{UInt32,9},GainMatrix4::Array{Float64,9},OldGainMatrix4::Array{Float64,9},GainTally4_K::AbstractArray{UInt32,9},OldGainTally4::Array{UInt32,9})

    # weighted average 
    @. OldGainMatrix3 = (GainMatrix3*GainTally3_K+OldGainMatrix3*OldGainTally3)/(GainTally3_K+OldGainTally3)
    @. OldGainMatrix4 = (GainMatrix4*GainTally4_K+OldGainMatrix4*OldGainTally4)/(GainTally4_K+OldGainTally4)

    replace!(OldGainMatrix3,NaN=>0e0)
    replace!(OldGainMatrix4,NaN=>0e0)

    # adding tallies 
    @. OldGainTally3 += GainTally3_K
    @. OldGainTally4 += GainTally4_K

end

"""
    WeightedAverageLossBinary!(LossMatrix,OldLossMatrix,LossTally,OldLossTally)

Computes the integral estimate by weighted average of the old and new gain matrices. Mutating the old gain and tally terms.
"""
function WeightedAverageLossBinary!(LossMatrix::Array{Float64,6},OldLossMatrix::Array{Float64,6},LossTally::Array{UInt32,6},OldLossTally::Array{UInt32,6})

    # weighted average 
    @. OldLossMatrix = (LossMatrix*LossTally+OldLossMatrix*OldLossTally)/(LossTally+OldLossTally)

    replace!(OldLossMatrix,NaN=>0e0)

    # adding tallies
    @. OldLossTally += LossTally

end