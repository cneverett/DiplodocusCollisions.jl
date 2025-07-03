
#= 
This module provides functions for MonteCarlo Integration of S and T Matrices
=#

"""
    BinaryMonteCarlo_Debug(Arrays,ArrayOfLocks,sigma,dsigmadt,Parameters,numLoss,numGain,scale,prog,thread_id)

Debug version of `BinaryMonteCarlo` which is only run when `numThreads=1`, useful for testing. Function performs the Monte-Carlo sampling of incoming and outgoing particle states for binary interactions.
"""
function BinaryMonteCarlo_Debug!(GainTotal3::Array{Float64,9},GainTotal4::Array{Float64,9},LossTotal::Array{Float64,6},GainTally3::Array{UInt32,9},GainTally4::Array{UInt32,9},LossTally::Array{UInt32,6},ArrayOfLocks,sigma::Function,dsigmadt::Function,Parameters::Tuple{String,String,String,String,Float64,Float64,Float64,Float64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64},numLoss::Int64,numGain::Int64,scale::Float64,prog::Progress,thread_id::Int64)

    # Set Parameters
    (name1,name2,name3,name4,m1,m2,m3,m4,p1_low,p1_up,p1_grid_st,p1_num,u1_grid_st,u1_num,h1_grid_st,h1_num,p2_low,p2_up,p2_grid_st,p2_num,u2_grid_st,u2_num,h2_grid_st,h2_num,p3_low,p3_up,p3_grid_st,p3_num,u3_grid_st,u3_num,h3_grid_st,h3_num,p4_low,p4_up,p4_grid_st,p4_num,u4_grid_st,u4_num,h4_grid_st,h4_num) = Parameters

    # Set up worker
    Threads.@spawn begin

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
    loc12::CartesianIndex{6} = CartesianIndex(0,0,0,0,0,0)

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

    LocalGainTotal3::Array{Float64,3} = zeros(Float64,size(GainTotal3)[1:3])
    LocalGainTally3::Array{UInt32,3} = zeros(UInt32,size(GainTally3)[1:3])
    LocalGainTotal4::Array{Float64,3} = zeros(Float64,size(GainTotal4)[1:3])
    LocalGainTally4::Array{UInt32,3} = zeros(UInt32,size(GainTally4)[1:3])

    @inbounds for _ in 1:numLoss
        
        # generate p1 and p2 vectors initially as to not have to re-calculate
        RPointSphereCosThetaPhi!(p1v)
        RPointSphereCosThetaPhi!(p2v)
        RPointLogMomentum!(p1v,p1_up,p1_low,p1_num)
        RPointLogMomentum!(p2v,p2_up,p2_low,p2_num)

        # LossVal
        (LossVal,sBig,sSmol) = LossValue(p1v,p2v,sigma,m1,m2,m3,m4)
        # Calculate T Array Location
        p1loc = location(p1_low,p1_up,p1_num,p1v[1],p1_grid)
        p2loc = location(p2_low,p2_up,p2_num,p2v[1],p2_grid)
        u1loc = location(u_low,u_up,u1_num,p1v[2],u1_grid)
        u2loc = location(u_low,u_up,u2_num,p2v[2],u2_grid)
        h1loc = location(h_low,h_up,h1_num,p1v[3],h1_grid)
        h2loc = location(h_low,h_up,h2_num,p2v[3],h2_grid)
        loc12 = CartesianIndex(p1loc,u1loc,h1loc,p2loc,u2loc,h2loc)

        fill!(LocalGainTally3,UInt32(0))
        fill!(LocalGainTally4,UInt32(0))

        if LossVal != 0e0 # i.e. it is a valid interaction state

            (w3,w4,t,h) = WeightedFactors(p1v,p2v,m1,m2,m3,m4,sBig,sSmol,scale)

            fill!(LocalGainTotal3,Float64(0))
            fill!(LocalGainTotal4,Float64(0))
                   
            @inbounds for _ in 1:numGain

                prob3 = RPointSphereWeighted!(p3v,w3)    
                prob4 = RPointSphereWeighted!(p4v,w4)
                RotateToLab!(p3v,p4v,t,h)
                @. p3pv = p3v
                @. p4pv = p4v
                
                # === p3 === #

                # Calculate p3 value
                (p_physical,pp_physical,NumStates) = MomentumValue!(p3v,p3pv,p1v,p2v,m1,m2,m3,m4,p3_low,p3_up)

                # Gain Array Tallies
                # For each u3,h3 sampled, p3 will be + or -ve, corresponding to a change in sign of u3 and a rotation of h3 by pi i.e. mod(h3+1,2). Therefore by sampling one u3,h3 we are actually sampling u3 and -u3 and h3, mod(h3+1,2) with one or both having valid p3 states. NOTE: This has been removed due to difference in sampling probability not being accounted for  
                #if NumStates != 0
                    u3loc = location(u_low,u_up,u3_num,p3v[2],u3_grid)
                    h3loc = location(h_low,h_up,h3_num,p3v[3],h3_grid)
                    LocalGainTally3[end,u3loc,h3loc] += UInt32(1)
                #end

                # Calculate Gain Array totals
                if NumStates == 1
                    if p_physical
                        p3loc = location(p3_low,p3_up,p3_num,p3v[1],p3_grid)
                        GainVal = GainValue3(p3v,p1v,p2v,sBig,sSmol,dsigmadt,m1,m2,m3,m4)
                        LocalGainTotal3[p3loc,u3loc,h3loc] += GainVal/prob3
                        LocalGainTally3[p3loc,u3loc,h3loc] += UInt32(1)
                    end
                end

                if NumStates == 2
                    if p_physical
                        p3loc = location(p3_low,p3_up,p3_num,p3v[1],p3_grid)
                        GainVal = GainValue3(p3v,p1v,p2v,sBig,sSmol,dsigmadt,m1,m2,m3,m4)
                        LocalGainTotal3[p3loc,u3loc,h3loc] += GainVal/prob3
                        LocalGainTally3[p3loc,u3loc,h3loc] += UInt32(1)
                    end
                    if pp_physical
                        u3ploc = location(u_low,u_up,u3_num,p3pv[2],u3_grid)
                        h3ploc = location(h_low,h_up,h3_num,p3pv[3],h3_grid)
                        p3ploc = location(p3_low,p3_up,p3_num,p3pv[1],p3_grid)
                        GainValp = GainValue3(p3pv,p1v,p2v,sBig,sSmol,dsigmadt,m1,m2,m3,m4)
                        LocalGainTotal3[p3ploc,u3ploc,h3ploc] += GainValp/prob3
                        LocalGainTally3[p3ploc,u3ploc,h3ploc] += UInt32(1)
                    end
                end

                # === p4 === #

                # Calculate p4 value
                (p_physical,pp_physical,NumStates) = MomentumValue!(p4v,p4pv,p2v,p1v,m2,m1,m4,m3,p4_low,p4_up)

                # S Array Tallies
                # For each u3,h4 sampled, p4 will be + or -ve, corresponding to a change in sign of u3 and a shift in h4 by pi i.e. Mod(h4+1,2). Therefore by sampling one u3 we are actually sampling u3/h4 and -u3/mod(h4+1,2) with one or both having valid p4 states. NOTE: This has been removed due to difference in sampling probability not being accounted for 
                u4loc = location(u_low,u_up,u4_num,p4v[2],u4_grid)
                h4loc = location(h_low,h_up,h4_num,p4v[3],h4_grid)
                #u4locMirror = location(u_low,u_up,u4_num,-p4v[2],u4_grid)
                #h4locMirror = location(h_low,h_up,h4_num,mod(p4v[3]+1e0,2e0),h4_grid)
                LocalGainTally4[end,u4loc,h4loc] += UInt32(1)
                #LocalGainTally4[end,u4locMirror,h4locMirror] += UInt32(1)

                # Calculate S Array totals
                if NumStates == 1
                    if p_physical
                        p4loc = location(p4_low,p4_up,p4_num,p4v[1],p4_grid)
                        GainVal = GainValue4(p4v,p1v,p2v,sBig,sSmol,dsigmadt,m1,m2,m3,m4,prob4)
                        LocalGainTotal4[p4loc,u4loc,h4loc] += GainVal/prob4
                        LocalGainTally4[p4loc,u4loc,h4loc] += UInt32(1)
                    end
                end

                if NumStates == 2
                    if p_physical
                        p4loc = location(p4_low,p4_up,p4_num,p4v[1],p4_grid)
                        GainVal = GainValue4(p4v,p1v,p2v,sBig,sSmol,dsigmadt,m1,m2,m3,m4,prob4)
                        LocalGainTotal4[p4loc,u4loc,h4loc] += GainVal/prob4
                        LocalGainTally4[p4loc,u4loc,h4loc] += UInt32(1)
                    end
                    if pp_physical
                        u4ploc = location(u_low,u_up,u4_num,p4pv[2],u4_grid)
                        h4ploc = location(h_low,h_up,h4_num,p4pv[3],h4_grid)
                        p4ploc = location(p4_low,p4_up,p4_num,p4pv[1],p4_grid)
                        GainValp = GainValue4(p4pv,p1v,p2v,sBig,sSmol,dsigmadt,m1,m2,m3,m4,prob4)
                        LocalGainTotal4[p4ploc,u4ploc,h4ploc] += GainValp/prob4
                        LocalGainTally4[p4ploc,u4ploc,h4ploc] += UInt32(1)
                    end
                end

            end # S loop

        else # no valid interaction state
            # add one to tally of all relevant S tallies i.e. all momenta and all angles as no emission states are possible
            @view(LocalGainTally3[end,:,:]) .+= UInt32(1)
            @view(LocalGainTally4[end,:,:]) .+= UInt32(1)
        end

        # assign values to arrays
        @lock ArrayOfLocks[p1loc] begin
            LossTotal[loc12] += LossVal
            LossTally[loc12] += UInt32(1)
            @view(GainTally3[:,:,:,loc12]) .+= LocalGainTally3
            @view(GainTally4[:,:,:,loc12]) .+= LocalGainTally4
            if LossVal != 0e0
                @view(GainTotal3[:,:,:,loc12]) .+= LocalGainTotal3
                @view(GainTotal4[:,:,:,loc12]) .+= LocalGainTotal4
            end
        end

        if thread_id == 1 # on main thread
            next!(prog)
        end

    end # T loop

    end # Thread spawn 

    return nothing

end # function 
