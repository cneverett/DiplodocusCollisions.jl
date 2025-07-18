"""
    EmissionInteractionIntegration(userSyncInputMultiThread)

Function to run the Monte Carlo integration of the S array in a serial environment. 
"""
function EmissionInteractionIntegration(Setup::Tuple{Tuple{String,String,String,String,Float64,Float64,Float64,Float64,Float64,Float64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Vector{Float64}},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64},Int64},Int64,Int64,Int64,String,String})

    # ======= Load User Parameters ======= #

    (Parameters,scale,numLoss,numGain,numThreads,fileLocation,fileName) = Setup;
    #(name1,name2,mu1,mu2,z1,z2,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,BMag) = Parameters;
    (name1,name2,name3,type,m1,m2,m3,z1,z2,z3,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num,Ext) = Parameters

    # ==================================== #

    # ========= Valid Grids? =============== #
        
        if u1_grid == "b" && u1_num%2 == 0 || u2_grid == "b" && u2_num%2 == 0
            error("Binary grid must have odd number of bins")
        end

    # ====================================== #

    # ======= Define Kernel Functions Based on Emission Type  ========= #

        name_kernel = Symbol(type*"Kernel")
        EmissionKernel = getfield(DiplodocusCollisions,name_kernel)

    # ============================================================================ #

    # ======== Set up Array of Locks ====== #

        ArrayOfLocks = [Threads.SpinLock() for _ in 1:p2_num]    

    # ===================================== #

    # ======== Load Old Arrays ========= #

        println("Loading Old and Allocating New Sampling Arrays")
        
        filePath = fileLocation*"\\"*fileName

        (OldGainTally2,OldGainTally3,OldLossTally1,OldGainMatrix2,OldGainMatrix3,OldLossMatrix1) = OldMonteCarloArraysEmission(Parameters,filePath)

        (GainTotal2,GainTotal3,LossTotal1,GainTally2,GainTally3,LossTally1,GainMatrix2,GainMatrix3,LossMatrix1) = MonteCarloArraysEmission(Parameters)

    # ================================= #

    # ===== Run MonteCarlo Integration ==== #

        println("Running Monte Carlo Integration")

        prog = Progress(numLoss)

        for (ii,scale_val) in enumerate(scale)
        numT = round(Int,numLoss/length(scale))

        println("")
        println("scale = $scale_val, itteration = $ii out of $(length(scale))")
        println("")

            # reset new arrays
            fill!(LossTotal1,Float64(0))
            fill!(GainTotal2,Float64(0))
            fill!(GainTotal3,Float64(0))
            fill!(LossTally1,UInt32(0))
            fill!(GainTally2,UInt32(0))
            fill!(GainTally3,UInt32(0))

            #workers  = [EmissionMonteCarloAxi_MultiThread!(SAtotal,SAtally,ArrayOfLocks,Parameters,numT,numSiterPerThread,nThreads,prog,thread) for thread in 1:nThreads]
            if numThreads == 1
                EmissionMonteCarlo_Debug!(GainTotal2,GainTally2,GainTotal3,GainTally3,LossTotal1,LossTally1,ArrayOfLocks,EmissionKernel,Parameters,numT,numGain,scale_val,prog,1)
            else 
                workers  = [EmissionMonteCarlo!(GainTotal2,GainTally2,GainTotal3,GainTally3,LossTotal1,LossTally1,ArrayOfLocks,EmissionKernel,Parameters,numT,numGain,scale_val,prog,thread) for thread in 1:numThreads]
                wait.(workers) # Allow all workers to finish
            end

    # ===================================== #


    # ===== Update Gain and Loss Matrices === #

            # Apply symmetries 
            GainLossSymmetryEmission!(GainTotal2,GainTotal3,GainTally2,GainTally3,LossTotal1,LossTally1)

            # Calculate Gain and Loss matrix
            GainMatrix2 = GainTotal2 ./ GainTally2
            GainMatrix3 = GainTotal3 ./ GainTally3
            LossMatrix1 = LossTotal1 ./ LossTally1

            println("")
            println("Applying Momentum Space Factors")

            # Angle / Momentum Ranges
            p3val = bounds(p3_low,p3_up,p3_num,p3_grid)
            u3val = bounds(u_low,u_up,u3_num,u3_grid)
            h3val = bounds(h_low,h_up,h3_num,h3_grid).*pi

            # Apply Momentum space volume elements
            MomentumSpaceFactorsEmission!(LossMatrix1,GainMatrix2,GainMatrix3,p3val,u3val,h3val)

            println("Weighting average of New and Old Sampling Arrays")

            # old arrays are modified in this process 
            WeightedAverageGainEmission!(GainMatrix2,OldGainMatrix2,GainTally2,OldGainTally2,GainMatrix3,OldGainMatrix3,GainTally3,OldGainTally3)
            WeightedAverageLossEmission!(LossMatrix1,OldLossMatrix1,LossTally1,OldLossTally1)

        end # scale loop 

        finish!(prog)

    # ===================================== #

    # ========== Save Arrays ============== #
            
#=         f = jldopen(filePath,"w") # creates file and overwrites previous file if one existed
        write(f,"STotal",SAtotal)
        write(f,"STally",SAtally)
        write(f,"SMatrix",SMatrix)
        #write(f,"pMax",pMax)
        #write(f,"tMinMax",tMinMax)
        write(f,"SConverge",SConverge)
        write(f,"Parameters",Parameters)
        close(f) =#

        println("")
        println("Saving Arrays")

        f = jldopen(filePath,"w") # creates file and overwrites previous file if one existed
        write(f,"GainTally2",OldGainTally2)
        write(f,"GainMatrix2",OldGainMatrix2)

        write(f,"GainTally3",OldGainTally3)
        write(f,"GainMatrix3",OldGainMatrix3)

        write(f,"LossTally1",OldLossTally1)
        write(f,"LossMatrix1",OldLossMatrix1)

        write(f,"Parameters",Parameters)

        close(f)

    # ===================================== #

end # function 