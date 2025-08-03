"""
    WeightedAverageGainBinary!(GainMatrix3,OldGainMatrix3,GainTally3_K,OldGainWeights3,GainMatrix4,OldGainMatrix4,GainTally4_K,OldGainWeights4_K)

Computes the integral estimate by weighted average of the old and new gain matrices. Mutating the old gain and tally terms.
```math
    I = (I1w1 + I2w2)/(w1 + w2)
```
where `I1` and `I1` are the old and new gain matrix element estimates and `w1` and `w2` are the corresponding weights. Here the weights are taken to be `w=k/N`
"""
function WeightedAverageGainBinary!(GainMatrix3::Array{Float64,9},OldGainMatrix3::Array{Float64,9},GainTally3_K::AbstractArray{UInt32,9},GainTally3_N::AbstractArray{UInt32,8},OldGainWeights3::Array{Float64,9},GainMatrix4::Array{Float64,9},OldGainMatrix4::Array{Float64,9},GainTally4_K::AbstractArray{UInt32,9},GainTally4_N::AbstractArray{UInt32,8},OldGainWeights4::Array{Float64,9})

    # new weights k/N
    NewGainWeights3 = similar(OldGainWeights3)
    NewGainWeights4 = similar(OldGainWeights4)
    for i in axes(GainTally3_K,1)
        @. NewGainWeights3[i,:,:,:,:,:,:,:,:] = GainTally3_K[i,:,:,:,:,:,:,:,:]^2 / GainTally3_N
    end
    for i in axes(GainTally4_K,1)
        @. NewGainWeights4[i,:,:,:,:,:,:,:,:] = GainTally4_K[i,:,:,:,:,:,:,:,:]^2 / GainTally4_N
    end
    replace!(NewGainWeights3,NaN=>0e0)
    replace!(NewGainWeights4,NaN=>0e0)

    # weighted average
    @. OldGainMatrix3 = (GainMatrix3*NewGainWeights3+OldGainMatrix3*OldGainWeights3)/(NewGainWeights3+OldGainWeights3)
    @. OldGainMatrix4 = (GainMatrix4*NewGainWeights4+OldGainMatrix4*OldGainWeights4)/(NewGainWeights4+OldGainWeights4)

    replace!(OldGainMatrix3,NaN=>0e0)
    replace!(OldGainMatrix4,NaN=>0e0)

    # adjusting new weights for next integration step
    @. OldGainWeights3 += NewGainWeights3
    @. OldGainWeights4 += NewGainWeights4

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