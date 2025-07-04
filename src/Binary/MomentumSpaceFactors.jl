"""
    GainLossSymmetryBinary(GainTotal3,GainTotal4,GainTally3,GainTally4,LossTotal,LossTally,Indistinguishable_34,m1,m2,m3,m4)

Applies various physical symmetries to the Gain and Loss terms for Binary (12->34) interactions to improve Monte Carlo sampling error. 
"""
function GainLossSymmetryBinary!(GainTotal3,GainTotal4,GainTally3,GainTally4,LossTotal,LossTally,m1,m2,m3,m4)

    # if the outgoing particles are indistinguishable or have identical masses (both conditions give m3==m4) then the Gain terms are identical assuming they have the same discretisation
    if m3 == m4
        @. GainTotal3 = GainTotal3 + GainTotal4
        @. GainTotal4 = GainTotal3
        @. GainTally3 = GainTally3 + GainTally4
        @. GainTally4 = GainTally3
    end

    # The Gain and Loss matrices are symmetric in two ways. 
    # FIRST: they are ALWAYS symmetric with respect to θ->π-θ for all particle momentum states
    # SECOND: if the incident masses are equal (m1==m2) then Gain and Loss are symmetric to swapping the incident particles
    
    tmp_total3 = zeros(Float64,size(GainTotal3)[1:2])
    tmp_total4 = zeros(Float64,size(GainTotal4)[1:2])
    tmp_tally3 = zeros(Float64,size(GainTally3)[1:2])
    tmp_tally4 = zeros(Float64,size(GainTally4)[1:2])

    # Particle 3 Gain terms
    for h2 in axes(GainTotal3,9), u2 in axes(GainTotal3,8), p2 in axes(GainTotal3,7), h1 in axes(GainTotal3,6), u1 in axes(GainTotal3,5), p1 in axes(GainTotal3,4), h3 in axes(GainTotal3,3)

        nu1 = size(GainTotal3)[5]
        nu2 = size(GainTotal3)[8]

        ViewGainTotal3 = @view(GainTotal3[:,1:end,h3,p1,u1,h1,p2,u2,h2])
        ViewGainTotal3Mirror123 = @view(GainTotal3[:,end:-1:1,h3,p1,nu1-u1+1,h1,p2,nu2-u2+1,h2])

        ViewGainTally3 = @view(GainTally3[:,1:end,h3,p1,u1,h1,p2,u2,h2])
        ViewGainTally3Mirror123 = @view(GainTally3[:,end:-1:1,h3,p1,nu1-u1+1,h1,p2,nu2-u2+1,h2])

        if m1 == m2 # Both first and second Symmetry true
            ViewGainTotal3Swap12 = @view(GainTotal3[:,1:end,h3,p2,u2,h2,p1,u1,h1])
            ViewGainTotal3Swap12Mirror123 = @view(GainTotal3[:,end:-1:1,h3,p2,nu2-u2+1,h2,p1,nu1-u1+1,h1])
            ViewGainTally3Swap12 = @view(GainTally3[:,1:end,h3,p2,u2,h2,p1,u1,h1])
            ViewGainTally3Swap12Mirror123 = @view(GainTally3[:,end:-1:1,h3,p2,nu2-u2+1,h2,p1,nu1-u1+1,h1])

            @. tmp_total3 = ViewGainTotal3 + ViewGainTotal3Mirror123 + ViewGainTotal3Swap12 + ViewGainTotal3Swap12Mirror123
            @. tmp_tally3 = ViewGainTally3 + ViewGainTally3Mirror123 + ViewGainTally3Swap12 + ViewGainTally3Swap12Mirror123

            @. ViewGainTotal3 = tmp_total3
            @. ViewGainTotal3Mirror123 = tmp_total3
            @. ViewGainTotal3Swap12 = tmp_total3
            @. ViewGainTotal3Swap12Mirror123 = tmp_total3

            @. ViewGainTally3 = tmp_tally3
            @. ViewGainTally3Mirror123 = tmp_tally3
            @. ViewGainTally3Swap12 = tmp_tally3
            @. ViewGainTally3Swap12Mirror123 = tmp_tally3

        else # only first symmetry true

            tmp_total3 = ViewGainTotal3 + ViewGainTotal3Mirror123
            tmp_tally3 = ViewGainTally3 + ViewGainTally3Mirror123

            @. ViewGainTotal3 = tmp_total3
            @. ViewGainTotal3Mirror123 = tmp_total3

            @. ViewGainTally3 = tmp_tally3
            @. ViewGainTally3Mirror123 = tmp_tally3

        end

    end

    # Particle 4 Gain terms
    for h2 in axes(GainTotal4,9), u2 in axes(GainTotal4,8), p2 in axes(GainTotal4,7), h1 in axes(GainTotal4,6), u1 in axes(GainTotal4,5), p1 in axes(GainTotal4,4), h4 in axes(GainTotal4,3)

        nu1 = size(GainTotal4)[5]
        nu2 = size(GainTotal4)[8]

        ViewGainTotal4 = @view(GainTotal4[:,1:end,h4,p1,u1,h1,p2,u2,h2])
        ViewGainTotal4Mirror123 = @view(GainTotal4[:,end:-1:1,h4,p1,nu1-u1+1,h1,p2,nu2-u2+1,h2])

        ViewGainTally4 = @view(GainTally4[:,1:end,h4,p1,u1,h1,p2,u2,h2])
        ViewGainTally4Mirror123 = @view(GainTally4[:,end:-1:1,h4,p1,nu1-u1+1,h1,p2,nu2-u2+1,h2])

        if m1 == m2 # Both first and second Symmetry true
            ViewGainTotal4Swap12 = @view(GainTotal4[:,1:end,h4,p2,u2,h2,p1,u1,h1])
            ViewGainTotal4Swap12Mirror123 = @view(GainTotal4[:,end:-1:1,h4,p2,nu2-u2+1,h2,p1,nu1-u1+1,h1])
            ViewGainTally4Swap12 = @view(GainTally4[:,1:end,h4,p2,u2,h2,p1,u1,h1])
            ViewGainTally4Swap12Mirror123 = @view(GainTally4[:,end:-1:1,h4,p2,nu2-u2+1,h2,p1,nu1-u1+1,h1])

            @. tmp_total4 = ViewGainTotal4 + ViewGainTotal4Mirror123 + ViewGainTotal4Swap12 + ViewGainTotal4Swap12Mirror123
            @. tmp_tally4 = ViewGainTally4 + ViewGainTally4Mirror123 + ViewGainTally4Swap12 + ViewGainTally4Swap12Mirror123

            @. ViewGainTotal4 = tmp_total4
            @. ViewGainTotal4Mirror123 = tmp_total4
            @. ViewGainTotal4Swap12 = tmp_total4
            @. ViewGainTotal4Swap12Mirror123 = tmp_total4

            @. ViewGainTally4 = tmp_tally4
            @. ViewGainTally4Mirror123 = tmp_tally4
            @. ViewGainTally4Swap12 = tmp_tally4
            @. ViewGainTally4Swap12Mirror123 = tmp_tally4

        else # only first symmetry true

            tmp_total4 = ViewGainTotal4 + ViewGainTotal4Mirror123
            tmp_tally4 = ViewGainTally4 + ViewGainTally4Mirror123

            @. ViewGainTotal4 = tmp_total4
            @. ViewGainTotal4Mirror123 = tmp_total4

            @. ViewGainTally4 = tmp_tally4
            @. ViewGainTally4Mirror123 = tmp_tally4

        end

    end

    # Particle Loss Terms
    for h2 in axes(LossTotal,6), u2 in axes(LossTotal,5), p2 in axes(LossTotal,4), h1 in axes(LossTotal,3), u1 in axes(LossTotal,2), p1 in axes(LossTotal,1)

        nu1 = size(LossTotal)[2]
        nu2 = size(LossTotal)[5]        

        if m1==m2 # Both first and second Symmetry true

            tmp_total = LossTotal[p1,u1,h1,p2,u2,h2] + LossTotal[p2,u2,h2,p1,u1,h1] + LossTotal[p1,nu1-u1+1,h1,p2,nu2-u2+1,h2] + LossTotal[p2,nu2-u2+1,h2,p1,nu1-u1+1,h1]
            tmp_tally = LossTally[p1,u1,h1,p2,u2,h2] + LossTally[p2,u2,h2,p1,u1,h1] + LossTally[p1,nu1-u1+1,h1,p2,nu2-u2+1,h2] + LossTally[p2,nu2-u2+1,h2,p1,nu1-u1+1,h1]

            LossTotal[p1,u1,h1,p2,u2,h2] = tmp_total
            LossTotal[p2,u2,h2,p1,u1,h1] = tmp_total
            LossTotal[p1,nu1-u1+1,h1,p2,nu2-u2+1,h2] = tmp_total
            LossTotal[p2,nu2-u2+1,h2,p1,nu1-u1+1,h1] = tmp_total

            LossTally[p1,u1,h1,p2,u2,h2] = tmp_tally
            LossTally[p2,u2,h2,p1,u1,h1] = tmp_tally
            LossTally[p1,nu1-u1+1,h1,p2,nu2-u2+1,h2] = tmp_tally
            LossTally[p2,nu2-u2+1,h2,p1,nu1-u1+1,h1] = tmp_tally

        else # only first symmetry true

            tmp_total = LossTotal[p1,u1,h1,p2,u2,h2] + LossTotal[p1,nu1-u1+1,h1,p2,nu2-u2+1,h2]
            tmp_tally = LossTally[p1,u1,h1,p2,u2,h2] + LossTally[p1,nu1-u1+1,h1,p2,nu2-u2+1,h2]

            LossTotal[p1,u1,h1,p2,u2,h2] = tmp_total
            LossTotal[p1,nu1-u1+1,h1,p2,nu2-u2+1,h2] = tmp_total

            LossTally[p1,u1,h1,p2,u2,h2] = tmp_tally
            LossTally[p1,nu1-u1+1,h1,p2,nu2-u2+1,h2] = tmp_tally

        end
    end 

end

"""
    MomentumSpaceFactorsBinary!(GainMatrix3,GainMatrix4,u3_bounds,u4_bounds,h3_bounds,h3_bounds,Indistinguishable_12)

Applies momentum space volume element to the Gain and Loss matrices such that they have the correct dimensions for use in particle transport, i.e. they have units of length^3/time (normalised by σT*c)
"""
function MomentumSpaceFactorsBinary!(GainMatrix3::Array{Float64,9},GainMatrix4::Array{Float64,9},u3_bounds::Vector{Float64},h3_bounds::Vector{Float64},u4_bounds::Vector{Float64},h4_bounds::Vector{Float64},Indistinguishable_12::Bool)

    # Momentum space volume elements
    for h2 in axes(GainMatrix3,9), u2 in axes(GainMatrix3,8), p2 in axes(GainMatrix3,7), h1 in axes(GainMatrix3,6), u1 in axes(GainMatrix3,5), p1 in axes(GainMatrix3,4) # common axes
        for h3 in axes(GainMatrix3,3), u3 in axes(GainMatrix3,2), p3 in 1:size(GainMatrix3,1)
            GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] *= (u3_bounds[u3+1]-u3_bounds[u3])*(h3_bounds[h3+1]-h3_bounds[h3]) # du3dh3
            GainMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2] /= (1e0+Float64(Indistinguishable_12))
        end
        for h4 in axes(GainMatrix4,3), u4 in axes(GainMatrix4,2), p4 in 1:size(GainMatrix4,1)
            GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] *= (u4_bounds[u4+1]-u4_bounds[u4])*(h4_bounds[h4+1]-h4_bounds[h4]) # du4dh4
            GainMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2] /= (1e0+Float64(Indistinguishable_12))
        end

    end

    return nothing

end
