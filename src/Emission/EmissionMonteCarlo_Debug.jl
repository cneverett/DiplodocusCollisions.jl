"""
    EmissionMonteCarloAxi_MultiThread!(SAtotal,SAtally,pMax,tMinMax,Parameters,numSiter)

# Arguments
- `GainTotal::Array{Float64,4}` : Array of stored integration totals for S matrix for 2+B->2+B+1 interaction
- `GainTally::Array{UInt32,4}` : Array of stored integration tallies for S matrix for 2+B->2+B+1 interaction
- `pMax::Array{Float64,3}` : Array of maximum momentum values for species 2
- `tMinMax::Array{Float64,3}` : Array of minimum and maximum theta values for species 2
- `Parameters::Tuple{Float64,Float64,Float64,Int64,Float64,Float64,Int64,Int64,Int64}` : Tuple of parameters for the interaction
- `numSiter::Int64` : Number of S iterations

# Output:
- Argument arrays GainTotal,GainTally are mutated to include the results of the Monte Carlo Integration.

# Calculation In Brief
- Random Sample points in each of these domains
    - RandomPointSphere for theta (for species 1,2)
    - RandomPointMomentum for p ( species 1,2)
- Take random points (p1,p2,t1,t2) and calculate Synchrotron emissivity
- Find position in S arrays and allocated tallies and totals accordingly.
"""
function EmissionMonteCarlo_Debug!(GainTotal2::Array{Float64,6},GainTally2::Array{UInt32,6},GainTotal3::Array{Float64,6},GainTally3::Array{UInt32,6},LossTotal1::Array{Float64,3},LossTally1::Array{UInt32,3},ArrayOfLocks,EmissionKernel::Function,Parameters::Tuple{String,String,String,String,Float64,Float64,Float64,Float64,Float64,Float64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Vector{Float64}},numLoss::Int64,numGain::Int64,scale::Float64,prog::Progress,thread_id::Int64)

    # Set Parameters
    (name1,name2,name3,type,m1,m2,m3,z1,z2,z3,p1_low,p1_up,p1_grid_st,p1_num,u1_grid_st,u1_num,h1_grid_st,h1_num,p2_low,p2_up,p2_grid_st,p2_num,u2_grid_st,u2_num,h2_grid_st,h2_num,p3_low,p3_up,p3_grid_st,p3_num,u3_grid_st,u3_num,h3_grid_st,h3_num,Ext) = Parameters

    # Set up workers
    #Threads.@spawn begin

    # allocate arrays
    p1v::Vector{Float64} = zeros(Float64,4)
    p2v::Vector{Float64} = zeros(Float64,4)
    p3v::Vector{Float64} = zeros(Float64,4)
    Sval::Float64 = 0e0
    prob::Float64 = 0e0
    #WeightFactors::Tuple{Float64,Float64,Float64} = (0e0,0e0,0e0)
    w::Float64 = 0e0
    t::Float64 = 0e0
    h::Float64 = 0e0

    p1_grid::GridType = Grid_String_to_Type(p1_grid_st)
    p2_grid::GridType = Grid_String_to_Type(p2_grid_st)
    p3_grid::GridType = Grid_String_to_Type(p3_grid_st)
    u1_grid::GridType = Grid_String_to_Type(u1_grid_st)
    u2_grid::GridType = Grid_String_to_Type(u2_grid_st)
    u3_grid::GridType = Grid_String_to_Type(u3_grid_st)
    h1_grid::GridType = Grid_String_to_Type(h1_grid_st)
    h2_grid::GridType = Grid_String_to_Type(h2_grid_st)
    h3_grid::GridType = Grid_String_to_Type(h3_grid_st)

    p1loc::Int64 = 0
    p2loc::Int64 = 0
    p3loc::Int64 = 0
    u1loc::Int64 = 0
    u2loc::Int64 = 0
    u3loc::Int64 = 0
    h1loc::Int64 = 0
    h2loc::Int64 = 0
    h3loc::Int64 = 0

    #SmallParameters = (p3_low,p3_up,p3_num,p3_grid,u3_num,u3_grid,h3_num,h3_grid,m1,m2,m3,z1,z2,z3,BMag)

    localGainTotal2::Array{Float64,3} = zeros(Float64,size(GainTotal2)[1:3])
    localGainTally2::Array{UInt32,3} = zeros(UInt32,size(GainTally2)[1:3])
    localGainTotal3::Array{Float64,3} = zeros(Float64,size(GainTotal3)[1:3])
    localGainTally3::Array{UInt32,3} = zeros(UInt32,size(GainTotal3)[1:3])

    for _ in 1:numLoss

        # generate p1v (emitting particle)
        RPointSphereCosThetaPhi!(p1v)
        RPointLogMomentum!(p1v,p1_low,p1_up,p1_num)
        p1loc = location(p1_low,p1_up,p1_num,p1v[1],p1_grid)
        u1loc = location(u_low,u_up,u1_num,p1v[2],u1_grid)
        h1loc = location(h_low,h_up,h1_num,p1v[3],h1_grid)

        #WeightFactors = WeightedFactorsEmission(p1v,m1,scale)
        (w,t,h) = WeightedFactorsEmission(p1v,m1,scale)

        fill!(localGainTotal3,Float64(0))
        fill!(localGainTally3,UInt32(0))

        for _ in 1:numGain

            #ImportanceSamplingSync!(p1v,p3v,localGainTally3,localGainTotal3,SmallParameters,WeightFactors)

            prob = RPointSphereWeighted!(p3v,w) # sample angles aligned to p1v
            RotateToLab!(p3v,t,h)   # rotate to z aligned
            RPointLogMomentum!(p3v,p3_low,p3_up,p3_num)

            # calculate S value
            Sval = EmissionKernel(p3v,p1v,m1,z1,Ext)

            # find S array location 
            p3loc = location(p3_low,p3_up,p3_num,p3v[1],p3_grid)
            u3loc = location(u_low,u_up,u3_num,p3v[2],u3_grid)
            h3loc = location(h_low,h_up,h3_num,p3v[3],h3_grid)

            localGainTally3[p3loc,u3loc,h3loc] += UInt32(1)
            localGainTotal3[p3loc,u3loc,h3loc] += Sval/prob

#=             # generate p1v (photon)
            RPointSphereCosTheta!(p1v)
            RPointLogMomentum!(p1v,p1_low,p1_up,p1_num)

            # calculate S value
            Sval = SyncKernel(p1v,p2v,m2,z2,BMag)
            # find S array location 
            p1loc = location(p1_low,p1_up,p1_num,p1v[1],p1_grid)
            u1loc = location(u_low,u_up,u1_num,p1v[2],u1_grid)

            localGainTally[p1loc,u1loc] += UInt32(1)
            localGainTotal[p1loc,u1loc] += Sval =#
        
        end

        # assign values to arrays
        @lock ArrayOfLocks[p1loc] begin
            @view(GainTally3[:,:,:,p1loc,u1loc,h1loc]) .+= localGainTally3
            @view(GainTotal3[:,:,:,p1loc,u1loc,h1loc]) .+= localGainTotal3
        end

        if thread_id == 1
            next!(prog)
        end

    end

    #end # workers 

end # function