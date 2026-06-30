"""
    BinaryMonteCarlo(Arrays,ArrayOfLocks,sigma,dsigmadt,Parameters,numLoss,numGain,scale,prog,thread_id)

Function performs the Monte-Carlo sampling of incoming and outgoing particle states for binary interactions.
"""
function BinaryMonteCarlo!(OldGainWeights3::ZArray,OldGainWeights4::ZArray,OldLossTally::ZArray,OldGainMatrix3::ZArray,OldGainMatrix4::ZArray,OldLossMatrix1::ZArray,OldLossMatrix2::ZArray,CorrectedGainMatrix3::ZArray,CorrectedGainMatrix4::ZArray,CorrectedLossMatrix1::ZArray,CorrectedLossMatrix2::ZArray,sigma::Function,dsigmadt::Function,Parameters::Tuple{String,String,String,String,Float64,Float64,Float64,Float64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64},numLoss::Int64,numGain::Int64,indices::Vector{CartesianIndex{2}},scale::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64},Int64},prog::Progress,thread_id::Int64)

    # Set Parameters
    (name1,name2,name3,name4,m1,m2,m3,m4,p1_low,p1_up,p1_grid_st,p1_num,u1_grid_st,u1_num,h1_grid_st,h1_num,p2_low,p2_up,p2_grid_st,p2_num,u2_grid_st,u2_num,h2_grid_st,h2_num,p3_low,p3_up,p3_grid_st,p3_num,u3_grid_st,u3_num,h3_grid_st,h3_num,p4_low,p4_up,p4_grid_st,p4_num,u4_grid_st,u4_num,h4_grid_st,h4_num) = Parameters

    # Set up worker
    Threads.@spawn begin

    println(stdout,"Thread ",thread_id," starting...")
    flush(stdout)

    # allocate arrays for each thread
    p1v::Vector{Float64} = zeros(Float64,4)
    p2v::Vector{Float64} = zeros(Float64,4)
    p3v::Vector{Float64} = zeros(Float64,4)
    p3pv::Vector{Float64} = zeros(Float64,4)
    p4v::Vector{Float64} = zeros(Float64,4)
    p4pv::Vector{Float64} = zeros(Float64,4)
    GainVal::Float64 = 0e0
    GainValp::Float64 = 0e0
    LossVal::Float64 = 0e0
    p_physical::Bool = true
    pp_physical::Bool = true
    NumStates::Int64 = 2
    prob3::Float64 = 0e0
    prob4::Float64 = 0e0
    w3::Float64 = 0e0
    w4::Float64 = 0e0
    t::Float64 = 0e0
    h::Float64 = 0e0

    p1loc::Int64 = 0
    p2loc::Int64 = 0
    u1loc::Int64 = 0
    u2loc::Int64 = 0
    h1loc::Int64 = 0
    h2loc::Int64 = 0
    p3loc::Int64 = 0
    u3loc::Int64 = 0
    h3loc::Int64 = 0
    p3ploc::Int64 = 0
    u3ploc::Int64 = 0
    h3ploc::Int64 = 0
    p4loc::Int64 = 0
    u4loc::Int64 = 0
    h4loc::Int64 = 0
    p4ploc::Int64 = 0
    u4ploc::Int64 = 0
    h4ploc::Int64 = 0
    loc12::CartesianIndex{4} = CartesianIndex(0,0,0,0)

    p1_grid::GridType = Grid_String_to_Type(p1_grid_st)
    p2_grid::GridType = Grid_String_to_Type(p2_grid_st)
    p3_grid::GridType = Grid_String_to_Type(p3_grid_st)
    p4_grid::GridType = Grid_String_to_Type(p4_grid_st)
    u1_grid::GridType = Grid_String_to_Type(u1_grid_st)
    u2_grid::GridType = Grid_String_to_Type(u2_grid_st)
    u3_grid::GridType = Grid_String_to_Type(u3_grid_st)
    u4_grid::GridType = Grid_String_to_Type(u4_grid_st)
    h1_grid::GridType = Grid_String_to_Type(h1_grid_st)
    h2_grid::GridType = Grid_String_to_Type(h2_grid_st)
    h3_grid::GridType = Grid_String_to_Type(h3_grid_st)
    h4_grid::GridType = Grid_String_to_Type(h4_grid_st)

    symmetric_grid::Bool = (u1_grid == "u" || u1_grid == "b") && (u2_grid == "u" || u2_grid == "b") && (u3_grid == "u" || u3_grid == "b") && (u4_grid == "u" || u4_grid == "b")

    Indistinguishable_12::Bool = name1 == name2
    Indistinguishable_34::Bool = name3 == name4

    u1_r::Vector{Float64} = bounds(u_low,u_up,u1_num,u1_grid_st)
    h1_r::Vector{Float64} = bounds(h_low,h_up,h1_num,h1_grid_st)
    u2_r::Vector{Float64} = bounds(u_low,u_up,u2_num,u2_grid_st)
    h2_r::Vector{Float64} = bounds(h_low,h_up,h2_num,h2_grid_st)

    u1_up::Float64 = 0.0
    u1_low::Float64 = 0.0
    h1_up::Float64 = 0.0
    h1_low::Float64 = 0.0
    u2_up::Float64 = 0.0
    u2_low::Float64 = 0.0
    h2_up::Float64 = 0.0
    h2_low::Float64 = 0.0

    # Angle / Momentum Ranges
    u3val = bounds(u_low,u_up,u3_num,u3_grid_st)
    u4val = bounds(u_low,u_up,u4_num,u4_grid_st)
    h3val = bounds(h_low,h_up,h3_num,h3_grid_st).*pi
    h4val = bounds(h_low,h_up,h4_num,h4_grid_st).*pi

    # local arrays are size of each chunk of Zarr that store data
    # for gain arrays this is p3_num+2,u3_num,h3_num,u1_num,h1_num,u2_num,h2_num
    # as p1 and p2 are what define the chunks
    # for loss arrays this is u1_num,h1_num,u2_num,h2_num as p1 and p2 are what define the chunks
    # these local arrays are in-memory and per thread
    # tallies for Gain terms are [k_underflow,k1,k2,k3,...,kn,k_overflow,N]
    ChunkGainTotal3::Array{Float64,7} = zeros(Float64,(p3_num+2),u3_num,h3_num,u1_num,h1_num,u2_num,h2_num)
    ChunkGainTally3::Array{UInt32,7} = zeros(UInt32,(p3_num+3),u3_num,h3_num,u1_num,h1_num,u2_num,h2_num)
    ChunkGainMatrix3::Array{Float64,7} = zeros(Float64,(p3_num+2),u3_num,h3_num,u1_num,h1_num,u2_num,h2_num)
    ChunkLossTotal::Array{Float64,4} = zeros(Float64,u1_num,h1_num,u2_num,h2_num)
    ChunkLossTally::Array{UInt32,4} = zeros(UInt32,u1_num,h1_num,u2_num,h2_num)
    if m3 != m4
        ChunkGainTotal4::Array{Float64,7} = zeros(Float64,(p4_num+2),u4_num,h4_num,u1_num,h1_num,u2_num,h2_num)
        ChunkGainTally4::Array{UInt32,7} = zeros(UInt32,(p4_num+3),u4_num,h4_num,u1_num,h1_num,u2_num,h2_num)
        ChunkGainMatrix4::Array{Float64,7} = zeros(Float64,(p4_num+2),u4_num,h4_num,u1_num,h1_num,u2_num,h2_num)
    end
    ChunkLossMatrix1::Array{Float64,4} = zeros(Float64,u1_num,h1_num,u2_num,h2_num)

    # local arrays for faster memory access that will be written to chunk arrays 
    LocalGainTotal3::Array{Float64,3} = zeros(Float64,size(ChunkGainTotal3)[1:3])
    LocalGainTally3::Array{UInt32,3} = zeros(UInt32,size(ChunkGainTally3)[1:3])
    if m3 != m4
        LocalGainTotal4::Array{Float64,3} = zeros(Float64,size(ChunkGainTotal4)[1:3])
        LocalGainTally4::Array{UInt32,3} = zeros(UInt32,size(ChunkGainTally4)[1:3])
    end

    # old chunk arrays 
    OldChunkGainMatrix3Full::Array{Float64,9} = zeros(Float64,(p3_num+2),u3_num,h3_num,1,u1_num,h1_num,1,u2_num,h2_num)
    OldChunkGainMatrix4Full::Array{Float64,9} = zeros(Float64,(p4_num+2),u4_num,h4_num,1,u1_num,h1_num,1,u2_num,h2_num)
    OldChunkLossMatrix1Full::Array{Float64,6} = zeros(Float64,1,u1_num,h1_num,1,u2_num,h2_num)
    OldChunkLossMatrix2Full::Array{Float64,6} = zeros(Float64,1,u2_num,h2_num,1,u1_num,h1_num)
    OldChunkGainWeights3Full::Array{Float64,9} = zeros(Float64,(p3_num+2),u3_num,h3_num,1,u1_num,h1_num,1,u2_num,h2_num)
    OldChunkGainWeights4Full::Array{Float64,9} = zeros(Float64,(p4_num+2),u4_num,h4_num,1,u1_num,h1_num,1,u2_num,h2_num)
    OldChunkLossTallyFull::Array{UInt32,6} = zeros(UInt32,1,u1_num,h1_num,1,u2_num,h2_num)

    OldChunkGainMatrix3 = @view(OldChunkGainMatrix3Full[:,:,:,1,:,:,1,:,:])
    OldChunkGainMatrix4 = @view(OldChunkGainMatrix4Full[:,:,:,1,:,:,1,:,:])
    OldChunkLossMatrix1 = @view(OldChunkLossMatrix1Full[1,:,:,1,:,:])
    OldChunkLossMatrix2 = @view(OldChunkLossMatrix2Full[1,:,:,1,:,:])
    OldChunkGainWeights3 = @view(OldChunkGainWeights3Full[:,:,:,1,:,:,1,:,:])
    OldChunkGainWeights4 = @view(OldChunkGainWeights4Full[:,:,:,1,:,:,1,:,:])
    OldChunkLossTally = @view(OldChunkLossTallyFull[1,:,:,1,:,:])

    for index in eachindex(indices)
        
        p1loc = indices[index][1]
        p2loc = indices[index][2]

        gain3loc = CartesianIndices((1:(p3_num+2), 1:u3_num, 1:h3_num,p1loc:p1loc, 1:u1_num, 1:h1_num,p2loc:p2loc, 1:u2_num, 1:h2_num))
        gain4loc = CartesianIndices((1:(p4_num+2), 1:u4_num, 1:h4_num,p1loc:p1loc, 1:u1_num, 1:h1_num,p2loc:p2loc, 1:u2_num, 1:h2_num))
        loss1loc = CartesianIndices((p1loc:p1loc, 1:u1_num, 1:h1_num,p2loc:p2loc, 1:u2_num, 1:h2_num))
        loss2loc = CartesianIndices((p2loc:p2loc, 1:u2_num, 1:h2_num,p1loc:p1loc, 1:u1_num, 1:h1_num))

        #println(stdout,"reading data on thread $thread_id for p1loc=$p1loc, p2loc=$p2loc")
        #flush(stdout)

        # Load old chunk arrays from Zarr
        #OldChunkGainMatrix3 .= OldGainMatrix3[:,:,:,p1loc,:,:,p2loc,:,:]
        #OldChunkGainMatrix4 .= OldGainMatrix4[:,:,:,p1loc,:,:,p2loc,:,:]
        #OldChunkLossMatrix1 .= OldLossMatrix1[p1loc,:,:,p2loc,:,:]
        #OldChunkGainWeights3 .= OldGainWeights3[:,:,:,p1loc,:,:,p2loc,:,:]
        #OldChunkGainWeights4 .= OldGainWeights4[:,:,:,p1loc,:,:,p2loc,:,:]
        #OldChunkLossTally .= OldLossTally[p1loc,:,:,p2loc,:,:]

        Zarr.readblock!(OldChunkGainMatrix3Full,OldGainMatrix3,gain3loc)
        Zarr.readblock!(OldChunkGainMatrix4Full,OldGainMatrix4,gain4loc)
        Zarr.readblock!(OldChunkLossMatrix1Full,OldLossMatrix1,loss1loc)
        #Zarr.readblock!(OldChunkLossMatrix2Full,OldLossMatrix2,loss2loc)
        Zarr.readblock!(OldChunkGainWeights3Full,OldGainWeights3,gain3loc)
        Zarr.readblock!(OldChunkGainWeights4Full,OldGainWeights4,gain4loc)
        Zarr.readblock!(OldChunkLossTallyFull,OldLossTally,loss1loc)

        #println(stdout,"read data on thread $thread_id for p1loc=$p1loc, p2loc=$p2loc")
        #flush(stdout)

        # reset in-memory local chunk arrays to zero 
        fill!(ChunkGainTotal3,Float64(0))
        fill!(ChunkGainTally3,UInt32(0))
        fill!(ChunkLossTotal,Float64(0))
        fill!(ChunkLossTally,UInt32(0))
        if m3 != m4
            fill!(ChunkGainTotal4,Float64(0))
            fill!(ChunkGainTally4,UInt32(0))
        end

        for (idx_scale,scale_val) in enumerate(scale)

        for u1loc in 1:u1_num, h1loc in 1:h1_num, u2loc in 1:u2_num, h2loc in 1:h2_num

            u1_up = u1_r[u1loc+1]
            u1_low = u1_r[u1loc]
            h1_up = h1_r[h1loc+1]
            h1_low = h1_r[h1loc]
            u2_up = u2_r[u2loc+1]
            u2_low = u2_r[u2loc]
            h2_up = h2_r[h2loc+1]
            h2_low = h2_r[h2loc]
            loc12 = CartesianIndex(u1loc,h1loc,u2loc,h2loc)

        for _ in 1:numLoss # sample incoming sates
        
            # generate p1 and p2 vectors initially as to not have to re-calculate
            RPointSphereCosThetaPhiBounds!(p1v,u1_low,u1_up,h1_low,h1_up)
            RPointSphereCosThetaPhiBounds!(p2v,u2_low,u2_up,h2_low,h2_up)
            RPointLogMomentum!(p1v,p1_up,p1_low,p1_num,p1loc)
            RPointLogMomentum!(p2v,p2_up,p2_low,p2_num,p2loc)

            fill!(LocalGainTally3,UInt32(0))
            if m3 != m4
                fill!(LocalGainTally4,UInt32(0))
            end

            # LossVal
            (LossVal,sBig,sSmol) = LossValue(p1v,p2v,sigma,m1,m2,m3,m4)

            if LossVal != 0e0 # i.e. it is a valid interaction state

                (w3,w4,t,h) = WeightedFactors(p1v,p2v,m1,m2,m3,m4,sBig,sSmol,scale_val)

                fill!(LocalGainTotal3,Float64(0))
                if m3 != m4
                    fill!(LocalGainTotal4,Float64(0))
                end
                    
                for _ in 1:(numGain*p3_num*u3_num*h3_num)

                    prob3 = RPointSphereWeighted!(p3v,w3)  
                    RotateToLab!(p3v,t,h)
                    @. p3pv = p3v

                    # Calculate p3 value
                    (p_physical,pp_physical,NumStates) = MomentumValue!(p3v,p3pv,p1v,p2v,m1,m2,m3,m4,p3_low,p3_up)

                    # Gain Array Tallies
                    # For each u3,h3 sampled, p3 will be + or -ve, corresponding to a change in sign of u3 and a rotation of h3 by pi i.e. mod(h3+1,2). Therefore by sampling one u3,h3 we are actually sampling u3 and -u3 and h3, mod(h3+1,2) with one or both having valid p3 states. NOTE: This has been removed due to difference in sampling probability not being accounted for  
                    u3loc = location(u_low,u_up,u3_num,p3v[2],u3_grid)
                    h3loc = location(h_low,h_up,h3_num,p3v[3],h3_grid)
                    LocalGainTally3[end,u3loc,h3loc] += UInt32(1)  # only need to do once even if there are two states 

                    # Calculate Gain Array totals
                    if NumStates == 1
                        if p_physical
                            p3loc = locationUnderOver(p3_low,p3_up,p3_num,p3v[1],p3_grid)
                            GainVal = GainValue3(p3v,p1v,p2v,sBig,sSmol,dsigmadt,m1,m2,m3,m4)
                            LocalGainTotal3[p3loc,u3loc,h3loc] += GainVal/prob3
                            LocalGainTally3[p3loc,u3loc,h3loc] += UInt32(1)
                        end
                    end

                    if NumStates == 2
                        if p_physical
                            p3loc = locationUnderOver(p3_low,p3_up,p3_num,p3v[1],p3_grid)
                            GainVal = GainValue3(p3v,p1v,p2v,sBig,sSmol,dsigmadt,m1,m2,m3,m4)
                            LocalGainTotal3[p3loc,u3loc,h3loc] += GainVal/prob3
                            LocalGainTally3[p3loc,u3loc,h3loc] += UInt32(1)
                        end
                        if pp_physical
                            u3ploc = location(u_low,u_up,u3_num,p3pv[2],u3_grid)
                            h3ploc = location(h_low,h_up,h3_num,p3pv[3],h3_grid)
                            p3ploc = locationUnderOver(p3_low,p3_up,p3_num,p3pv[1],p3_grid)
                            GainValp = GainValue3(p3pv,p1v,p2v,sBig,sSmol,dsigmadt,m1,m2,m3,m4)
                            LocalGainTotal3[p3ploc,u3ploc,h3ploc] += GainValp/prob3
                            LocalGainTally3[p3ploc,u3ploc,h3ploc] += UInt32(1)
                            #LocalGainTally3[end,u3ploc,h3ploc] += UInt32(1)
                        end
                    end

                end # p3 outgoing loop

                if m3 != m4

                    for _ in 1:(numGain*p4_num*u4_num*h4_num)

                        prob4 = RPointSphereWeighted!(p4v,w4)
                        RotateToLab!(p4v,t,h)
                        @. p4pv = p4v

                        # Calculate p4 value
                        (p_physical,pp_physical,NumStates) = MomentumValue!(p4v,p4pv,p2v,p1v,m2,m1,m4,m3,p4_low,p4_up)

                        # S Array Tallies
                        # For each u3,h4 sampled, p4 will be + or -ve, corresponding to a change in sign of u3 and a shift in h4 by pi i.e. Mod(h4+1,2). Therefore by sampling one u3 we are actually sampling u3/h4 and -u3/mod(h4+1,2) with one or both having valid p4 states. NOTE: This has been removed due to difference in sampling probability not being accounted for 
                        u4loc = location(u_low,u_up,u4_num,p4v[2],u4_grid)
                        h4loc = location(h_low,h_up,h4_num,p4v[3],h4_grid)
                        LocalGainTally4[end,u4loc,h4loc] += UInt32(1) # only need to do once even if there are two states 

                        # Calculate S Array totals
                        if NumStates == 1
                            if p_physical
                                p4loc = locationUnderOver(p4_low,p4_up,p4_num,p4v[1],p4_grid)
                                GainVal = GainValue4(p4v,p1v,p2v,sBig,sSmol,dsigmadt,m1,m2,m3,m4)
                                LocalGainTotal4[p4loc,u4loc,h4loc] += GainVal/prob4
                                LocalGainTally4[p4loc,u4loc,h4loc] += UInt32(1)
                            end
                        end

                        if NumStates == 2
                            if p_physical
                                p4loc = locationUnderOver(p4_low,p4_up,p4_num,p4v[1],p4_grid)
                                GainVal = GainValue4(p4v,p1v,p2v,sBig,sSmol,dsigmadt,m1,m2,m3,m4)
                                LocalGainTotal4[p4loc,u4loc,h4loc] += GainVal/prob4
                                LocalGainTally4[p4loc,u4loc,h4loc] += UInt32(1)
                            end
                            if pp_physical
                                u4ploc = location(u_low,u_up,u4_num,p4pv[2],u4_grid)
                                h4ploc = location(h_low,h_up,h4_num,p4pv[3],h4_grid)
                                p4ploc = locationUnderOver(p4_low,p4_up,p4_num,p4pv[1],p4_grid)
                                GainValp = GainValue4(p4pv,p1v,p2v,sBig,sSmol,dsigmadt,m1,m2,m3,m4)
                                LocalGainTotal4[p4ploc,u4ploc,h4ploc] += GainValp/prob4
                                LocalGainTally4[p4ploc,u4ploc,h4ploc] += UInt32(1)
                            end
                        end

                    end # p4 outgoing loop
            
                end

            else # no valid interaction state
                # add one to tally of all relevant S tallies i.e. all momenta and all angles as no emission states are possible
                @view(LocalGainTally3[end,:,:]) .+= UInt32(1)
                if m3 != m4
                    @view(LocalGainTally4[end,:,:]) .+= UInt32(1)
                end
            end

            # assign local arrays to chunk arrays
            ChunkLossTotal[loc12] += LossVal
            ChunkLossTally[loc12] += UInt32(1)
            @view(ChunkGainTally3[:,:,:,loc12]) .+= LocalGainTally3
            if m3 != m4 
                @view(ChunkGainTally4[:,:,:,loc12]) .+= LocalGainTally4
            end
            if LossVal != 0e0
                @view(ChunkGainTotal3[:,:,:,loc12]) .+= LocalGainTotal3
                if m3 != m4
                    @view(ChunkGainTotal4[:,:,:,loc12]) .+= LocalGainTotal4
                end
            end

        end # numLoss loop

        end # u1,h1,u2,h2 loop

        # === Update Gain and Loss Matrices === #

            # N values are last element of the tally array
            ChunkGainTally3_N = @view(ChunkGainTally3[end,:,:,:,:,:,:])
            # K value are all but last element of the tally array
            ChunkGainTally3_K = @view(ChunkGainTally3[1:end-1,:,:,:,:,:,:])
            if m3 != m4 
                ChunkGainTally4_N = @view(ChunkGainTally4[end,:,:,:,:,:,:])
                ChunkGainTally4_K = @view(ChunkGainTally4[1:end-1,:,:,:,:,:,:])
            end 

            # calculate the gain and loss matrices
            for i in axes(ChunkGainTotal3,1)
                @view(ChunkGainMatrix3[i,:,:,:,:,:,:]) .= @view(ChunkGainTotal3[i,:,:,:,:,:,:]) ./ ChunkGainTally3_N
            end
            replace!(ChunkGainMatrix3,NaN=>0e0); # remove NaN caused by / 0
            if m3 != m4
                for i in axes(ChunkGainTotal4,1)
                    @view(ChunkGainMatrix4[i,:,:,:,:,:,:]) .= @view(ChunkGainTotal4[i,:,:,:,:,:,:]) ./ ChunkGainTally4_N
                end
                replace!(ChunkGainMatrix4,NaN=>0e0); # remove NaN caused by /0e0
            end
            @. ChunkLossMatrix1 = ChunkLossTotal / ChunkLossTally;
            replace!(ChunkLossMatrix1,NaN=>0e0);

            # Momentum space volume elements
            if m3 == m4
                MomentumSpaceFactorsBinaryChunk!(ChunkGainMatrix3,u3val,h3val,Indistinguishable_12)
            else
                MomentumSpaceFactorsBinaryChunk!(ChunkGainMatrix3,ChunkGainMatrix4,u3val,h3val,u4val,h4val,Indistinguishable_12)
            end

            # perform weighted average of old and new gain matrices and loss matrices
            # old arrays are modified in this process
            if m3 == m4
                WeightedAverageGainBinaryChunk!(ChunkGainMatrix3,OldChunkGainMatrix3,ChunkGainTally3_K,ChunkGainTally3_N,OldChunkGainWeights3)
                if Indistinguishable_34 == false # particles are distinguishable 
                    @. OldChunkGainMatrix4 = OldChunkGainMatrix3
                    @. OldChunkGainWeights4 = OldChunkGainWeights3
                end
            else
                WeightedAverageGainBinaryChunk!(ChunkGainMatrix3,OldChunkGainMatrix3,ChunkGainTally3_K,ChunkGainTally3_N,OldChunkGainWeights3,ChunkGainMatrix4,OldChunkGainMatrix4,ChunkGainTally4_K,ChunkGainTally4_N,OldChunkGainWeights4)
            end
            WeightedAverageLossBinaryChunk!(ChunkLossMatrix1,OldChunkLossMatrix1,ChunkLossTally,OldChunkLossTally)


        end # scale loop  
        
        # ========= Apply Symmetries ========== # 

            # Apply Symmetries to the Gain and Loss Matrices, this does not affect the weighting of the average, which is already done in the previous step
            # TODO: improve this averaging to account for the weights and tallies to do a weighted average symmetrising
            GainLossPolarSymmetryMatrixBinaryChunk!(OldChunkGainMatrix3,OldChunkGainMatrix4,OldChunkLossMatrix1,OldChunkGainWeights3,OldChunkGainWeights4,OldChunkLossTally,m1,m2,m3,m4,symmetric_grid)
            GainLossAzimuthalSymmetryMatrixBinaryChunk!(OldChunkGainMatrix3,OldChunkGainMatrix4,OldChunkLossMatrix1,OldChunkGainWeights3,OldChunkGainWeights4,OldChunkLossTally,m1,m2,m3,m4)
            replace!(OldChunkGainMatrix3,NaN=>0e0); # remove NaN caused by /0e0 if weights are zero
            replace!(OldChunkGainMatrix4,NaN=>0e0); # remove NaN caused by /0e0 if weights are zero
            replace!(OldChunkLossMatrix1,NaN=>0e0); # remove NaN caused by /0e0 if tallies are zero

        # =========== Generate LossMatrix2 ==== #

            if Indistinguishable_12 == false # particles are distinguishable
                perm = [3,4,1,2]
                OldChunkLossMatrix2 .= permutedims(OldChunkLossMatrix1,perm)
            else
                fill!(OldChunkLossMatrix2,Float64(0))
            end

        # ===== Saving Uncorrected Arrays ===== #



            Zarr.writeblock!(OldChunkGainMatrix3Full,OldGainMatrix3,gain3loc)
            Zarr.writeblock!(OldChunkGainMatrix4Full,OldGainMatrix4,gain4loc)
            Zarr.writeblock!(OldChunkLossMatrix1Full,OldLossMatrix1,loss1loc)
            Zarr.writeblock!(OldChunkLossMatrix2Full,OldLossMatrix2,loss2loc)
            Zarr.writeblock!(OldChunkGainWeights3Full,OldGainWeights3,gain3loc)
            Zarr.writeblock!(OldChunkGainWeights4Full,OldGainWeights4,gain4loc)
            Zarr.writeblock!(OldChunkLossTallyFull,OldLossTally,loss1loc)

            #println(stdout,"written uncorrected data on thread $thread_id for p1loc=$p1loc, p2loc=$p2loc")
            #flush(stdout)

        # ===== Generate Corrected Arrays ===== #

            # modifies uncorrected arrays to produce corrected arrays, which are then saved to Zarr
            GainCorrectionChunk!(Parameters,OldChunkGainMatrix3,OldChunkGainMatrix4,OldChunkLossMatrix1,OldChunkLossMatrix2,p1loc,p2loc)

        # ========== Save Chunks to Zarr ============== #

            

            #OldGainMatrix3[:,:,:,p1loc,:,:,p2loc,:,:] = OldChunkGainMatrix3
            #OldGainMatrix4[:,:,:,p1loc,:,:,p2loc,:,:] = OldChunkGainMatrix4
            #OldLossMatrix1[p1loc,:,:,p2loc,:,:] = OldChunkLossMatrix1
            #OldLossMatrix2[p2loc,:,:,p1loc,:,:] = OldChunkLossMatrix2
            #OldGainWeights3[:,:,:,p1loc,:,:,p2loc,:,:] = OldChunkGainWeights3
            #OldGainWeights4[:,:,:,p1loc,:,:,p2loc,:,:] = OldChunkGainWeights4
            #OldLossTally[p1loc,:,:,p2loc,:,:] = OldChunkLossTally

            #CorrectedGainMatrix3[:,:,:,p1loc,:,:,p2loc,:,:] = CorrectedChunkGainMatrix3
            #CorrectedGainMatrix4[:,:,:,p1loc,:,:,p2loc,:,:] = CorrectedChunkGainMatrix4
            #CorrectedLossMatrix1[p1loc,:,:,p2loc,:,:] = CorrectedChunkLossMatrix1
            #CorrectedLossMatrix2[p2loc,:,:,p1loc,:,:] = CorrectedChunkLossMatrix2

            #println(stdout,"writing corrected data on thread $thread_id for p1loc=$p1loc, p2loc=$p2loc")
            #flush(stdout)

            Zarr.writeblock!(OldChunkGainMatrix3Full,CorrectedGainMatrix3,gain3loc)
            Zarr.writeblock!(OldChunkGainMatrix4Full,CorrectedGainMatrix4,gain4loc)
            Zarr.writeblock!(OldChunkLossMatrix1Full,CorrectedLossMatrix1,loss1loc)
            Zarr.writeblock!(OldChunkLossMatrix2Full,CorrectedLossMatrix2,loss2loc)

            println(stdout,"Completed MC loop on thread $thread_id for p1loc=$p1loc, p2loc=$p2loc")
            flush(stdout)

            # Update progress 
            next!(prog)

    end # indices loop

        println(stdout,"Thread ",thread_id," finished")
        flush(stdout)

    end # Thread spawn 

end # function 
