#= Script for running the ST integration and returning data arrays =#

#=
    Optimisation of the multi-threaded version would not have been possible without the kind guidance of those on the Julia Forums: https://discourse.julialang.org/t/fast-multi-threaded-array-editing-without-data-races/114863/41
    In particular the assistance of users: mbauman, adienes, Oscar_Smith, Satvik, Salmon, sgaure and foobar_lv2
=#

"""
    BinaryInteractionIntegration(Setup)

Function to build/load the collision matrices, then run the Monte Carlo integration of the Gain and Loss arrays in a multi-threaded environment. The function will run the Monte Carlo integration in parallel across the number of threads specified in the global variable `numThreads``. The function will then evaluate the Gain and Loss matrices using a weighted averaging method and save the results to a file specified by `fileLocation` and `fileName` in the input `Setup`.
"""
function BinaryInteractionIntegration(Setup::Tuple{Tuple{String,String,String,String,Float64,Float64,Float64,Float64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64},Tuple{Int64,Int64,Int64,Int64},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64},Int64},Int64,Int64,Int64,String,String})

    # ========= Load user Parameters ======= #

    (Parameters,bins,scale,numLoss,numGain,numThreads,fileLocation,fileName) = Setup
    (name1,name2,name3,name4,m1,m2,m3,m4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num,h4_grid,h4_num) = Parameters
    (p1loc_low,p1loc_up,p2loc_low,p2loc_up) = bins 

    # ====================================== #

    # ========= Valid Grids? =============== #
        
        if u1_grid == "b" && u1_num%2 == 0 || u2_grid == "b" && u2_num%2 == 0 || u3_grid == "b" && u3_num%2 == 0 || u4_grid == "b" && u4_num%2 == 0
            error("Binary grid must have odd number of bins")
        end

    # ====================================== #

    # ===== Are states Distinguishable ===== #

        Indistinguishable_12::Bool = name1 == name2
        Indistinguishable_34::Bool = name3 == name4

    # ====================================== #

    # ======= Define Cross Section Functions Based on Particle Selections ========= #

        name_sigma = Symbol("sigma_"*name1*name2*name3*name4)
        sigma = getfield(DiplodocusCollisions,name_sigma)
        name_dsigmadt = Symbol("dsigmadt_"*name1*name2*name3*name4)
        dsigmadt = getfield(DiplodocusCollisions,name_dsigmadt)

    # ============================================================================ #

    # ======== Set up Array of Locks ====== #

        ArrayOfLocks = [Threads.SpinLock() for _ in 1:p1_num] 

    # ===================================== #

    # =========== Load Old and Generate New Arrays ========= #

        println("Loading Old and Allocating New Sampling Arrays")
                
        filePath = fileLocation*"\\"*fileName

        (OldGainWeights3,OldGainWeights4,OldLossTally,OldGainMatrix3,OldGainMatrix4,OldLossMatrix1,OldLossMatrix2) = OldMonteCarloArraysBinary(Parameters,filePath)

        (GainTotal3,GainTotal4,LossTotal,GainTally3,GainTally4,LossTally,GainMatrix3,GainMatrix4,LossMatrix1,LossMatrix2) = MonteCarloArraysBinary(Parameters)

    # ===================================== #

    # ===== Run MonteCarlo Integration ==== #

        println("Running Monte Carlo Integration")

        for (ii,scale_val) in enumerate(scale)

            println("")
            println("scale = $scale_val, itteration = $ii out of $(length(scale))")
            println("")

            indices = CartesianIndices((p1loc_low:p1loc_up,p2loc_low:p2loc_up))
            length_indices::Int64 = length(indices)

            if length_indices/numThreads > 1.0
                length_div_threads::Int64 = ceil(Int64,length_indices/numThreads)
                index_range = Vector(0:length_div_threads:numThreads*length_div_threads)
                index_range[end] = length_indices
            else 
                index_range = 0:length_indices
            end

            # reset arrays
            fill!(GainTotal3,Float64(0))
            fill!(LossTotal,Float64(0))
            fill!(GainTally3,UInt32(0))
            fill!(LossTally,UInt32(0))

            fill!(GainTotal4,Float64(0))
            fill!(GainTally4,UInt32(0))

            #if numThreads == 1
            #    # Run in serial if only one thread, easier to use for debugging
            #    BinaryMonteCarlo_Debug!(GainTotal3,GainTotal4,LossTotal,GainTally3,GainTally4,LossTally,ArrayOfLocks,sigma,dsigmadt,Parameters,numT,numGain,scale_val,prog,1)
            #else
            #    workers = [BinaryMonteCarlo!(GainTotal3,GainTotal4,LossTotal,GainTally3,GainTally4,LossTally,ArrayOfLocks,sigma,dsigmadt,Parameters,numT,numGain,scale_val,prog,thread) for thread in 1:numThreads]
            #    wait.(workers) # Allow all workers to finish
            #end

            if numThreads == 1
                numProgress = numLoss*index_range[end]*u1_num*h1_num*u2_num*h2_num
                prog = Progress(numProgress)
                # Run in serial if only one thread, easier to use for debugging
                BinaryMonteCarlo_Debug!(GainTotal3,GainTotal4,LossTotal,GainTally3,GainTally4,LossTally,ArrayOfLocks,sigma,dsigmadt,Parameters,numLoss,numGain,indices,scale_val,prog,1)
            else
                numProgress = numLoss*index_range[1+1]*u1_num*h1_num*u2_num*h2_num
                prog = Progress(numProgress)
                workers = [BinaryMonteCarlo!(GainTotal3,GainTotal4,LossTotal,GainTally3,GainTally4,LossTally,ArrayOfLocks,sigma,dsigmadt,Parameters,numLoss,numGain,indices[index_range[thread]+1:index_range[thread+1]],scale_val,prog,thread) for thread in 1:(length(index_range)-1)]
                wait.(workers) # Allow all workers to finish
            end

            #thread =1

            #BinaryMonteCarlo!(GainTotal3,GainTotal4,LossTotal,GainTally3,GainTally4,LossTally,ArrayOfLocks,sigma,dsigmadt,Parameters,numLoss,numGain,bins,scale_val,thread)

    
    # ===================================== #

    # === Update Gain and Loss Matrices === #

            # N values are last element of the tally array
            GainTally3_N = @view(GainTally3[end,:,:,:,:,:,:,:,:])
            # K value are all but last element of the tally array
            GainTally3_K = @view(GainTally3[1:end-1,:,:,:,:,:,:,:,:])
            if m3 != m4 
                GainTally4_N = @view(GainTally4[end,:,:,:,:,:,:,:,:])
                GainTally4_K = @view(GainTally4[1:end-1,:,:,:,:,:,:,:,:])
            end

            println("")
            println("Applying Symmetries")

            # Apply Symmetries to the Gain and Loss Totals and Tallies
            GainLossSymmetryBinary!(GainTotal3,GainTotal4,GainTally3,GainTally4,LossTotal,LossTally,m1,m2,m3,m4)

            println("Generating New Sampling Arrays")

            # calculate the gain and loss matrices
            if m3 == m4
                # Only need to calculate GainMatrix3
                for i in axes(GainTotal3,1)
                    @view(GainMatrix3[i,:,:,:,:,:,:,:,:]) .= @view(GainTotal3[i,:,:,:,:,:,:,:,:]) ./ GainTally3_N
                end
                replace!(GainMatrix3,NaN=>0e0); # remove NaN caused by / 0
                # GainMatrix4 remains zero
            else
                for i in axes(GainTotal3,1)
                    @view(GainMatrix3[i,:,:,:,:,:,:,:,:]) .= @view(GainTotal3[i,:,:,:,:,:,:,:,:]) ./ GainTally3_N
                end
                replace!(GainMatrix3,NaN=>0e0); # remove NaN caused by /0e0
                for i in axes(GainTotal4,1)
                    @view(GainMatrix4[i,:,:,:,:,:,:,:,:]) .= @view(GainTotal4[i,:,:,:,:,:,:,:,:]) ./ GainTally4_N
                end
                replace!(GainMatrix4,NaN=>0e0); # remove NaN caused by /0e0
            end
            @. LossMatrix1 = LossTotal / LossTally;
            replace!(LossMatrix1,NaN=>0e0);

            println("Applying Momentum Space Factors")

            # Angle / Momentum Ranges
            u3val = bounds(u_low,u_up,u3_num,u3_grid)
            u4val = bounds(u_low,u_up,u4_num,u4_grid)
            h3val = bounds(h_low,h_up,h3_num,h3_grid).*pi
            h4val = bounds(h_low,h_up,h4_num,h4_grid).*pi

            # Momentum space volume elements
            if m3 == m4
                MomentumSpaceFactorsBinary!(GainMatrix3,u3val,h3val,Indistinguishable_12)
            else
                MomentumSpaceFactorsBinary!(GainMatrix3,GainMatrix4,u3val,h3val,u4val,h4val,Indistinguishable_12)
            end
                                        
            println("Weighting average of New and Old Sampling Arrays")

            # old arrays are modified in this process
            if m3 == m4
                WeightedAverageGainBinary!(GainMatrix3,OldGainMatrix3,GainTally3_K,GainTally3_N,OldGainWeights3)
                if Indistinguishable_34 == false # particles are distinguishable 
                    @. OldGainMatrix4 = OldGainMatrix3
                    @. OldGainWeights4 = OldGainWeights3
                end
            else
                WeightedAverageGainBinary!(GainMatrix3,OldGainMatrix3,GainTally3_K,GainTally3_N,OldGainWeights3,GainMatrix4,OldGainMatrix4,GainTally4_K,GainTally4_N,OldGainWeights4)
            end
            WeightedAverageLossBinary!(LossMatrix1,OldLossMatrix1,LossTally,OldLossTally)

        end # scale loop 

        if Indistinguishable_12 == false # particles are distinguishable
            perm = [4,5,6,1,2,3]
            OldLossMatrix2 .= permutedims(OldLossMatrix1,perm)
        else
            fill!(OldLossMatrix2,Float64(0))
        end

    # ===================================== #

    # ============= Error Estimates ======= # 

        println("Calculating Error Estimates")

        ErrorOutput =  DoesConserve((Parameters,OldGainMatrix3,OldGainMatrix4,OldLossMatrix1,OldLossMatrix2))

    # ===================================== #

    # ============= Error Estimates ======= # 

        println("Calculating Noise Corrected Arrays")

        CorrectedGainMatrix3, CorrectedGainMatrix4, CorrectedLossMatrix1, CorrectedLossMatrix2 =  GainCorrection(Parameters,OldGainMatrix3,OldGainMatrix4,OldLossMatrix1,OldLossMatrix2)

    # ===================================== #

    # =========== Garbage Collection ====== #

        (GainTotal3,GainTotal4,LossTotal,GainTally3,GainTally4,LossTally,GainMatrix3,GainMatrix4,LossMatrix1,LossMatrix2) = (0,0,0,0,0,0,0,0,0,0)
        GC.gc()
        
    # ===================================== #

    # ========== Save Arrays ============== #

        println("Saving Arrays")

        f = jldopen(filePath,"w") # creates file and overwrites previous file if one existed
        write(f,"GainWeights3",OldGainWeights3)
        write(f,"GainMatrix3",OldGainMatrix3)

        write(f,"GainWeights4",OldGainWeights4)
        write(f,"GainMatrix4",OldGainMatrix4)

        write(f,"LossTally",OldLossTally)
        write(f,"LossMatrix1",OldLossMatrix1)
        write(f,"LossMatrix2",OldLossMatrix2)

        write(f,"Parameters",Parameters)
        write(f,"ErrorEstimates",ErrorOutput)

        write(f,"CorrectedGainMatrix3",CorrectedGainMatrix3)
        write(f,"CorrectedGainMatrix4",CorrectedGainMatrix4)
        write(f,"CorrectedLossMatrix1",CorrectedLossMatrix1)
        write(f,"CorrectedLossMatrix2",CorrectedLossMatrix2)
        close(f)

    # ===================================== #

    GC.gc()

    return nothing

end # function