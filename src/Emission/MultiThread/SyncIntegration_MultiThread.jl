"""
    SpectraEvaluateMultiThreadEmission(userSyncInputMultiThread)

Function to run the Monte Carlo integration of the S array in a serial environment. 
"""
function SpectraEvaluateMultiThreadEmission(userInputEmissionMultiThread::Tuple{Tuple{String,String,String,String,Float64,Float64,Float64,Float64,Float64,Float64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64},Int64,Int64,Int64,String,String})

    # ======= Load User Parameters ======= #

    (Parameters,numTiterPerThread,numSiterPerThread,nThreads,fileLocation,fileName) = userInputEmissionMultiThread;
    #(name1,name2,mu1,mu2,z1,z2,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,BMag) = Parameters;
    (name1,name2,name3,type,mu1,mu2,mu3,z1,z2,z3,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num,BMag) = Parameters

    # ==================================== #

    # ========= Valid Grids? =============== #
        
        if u1_grid == "b" && u1_num%2 == 0 || u2_grid == "b" && u2_num%2 == 0
            error("Binary grid must have odd number of bins")
        end

    # ====================================== #

    # ======== Set up Array of Locks ====== #

        ArrayOfLocks = [Threads.SpinLock() for _ in 1:p2_num]    

    # ===================================== #

    # ======== Load Old Arrays ========= #

        println("Loading Old and Allocating New Sampling Arrays")
        
        filePath = fileLocation*"\\"*fileName

        (OldGainTally2,OldGainTally3,OldLossTally1,OldGainMatrix2,OldGainMatrix3,OldLossMatrix1) = OldMonteCarloArraysEmission(Parameters,filePath)

        (GainTotal2,GainTotal3,LossTotal1,GainTally2,GainTally3,LossTally1,GainMatrix2,GainMatrix3,LossMatrix1) = MonteCarloArraysEmission(Parameters)

#=         if fileExist
            f = jldopen(filePath,"r+");
            SAtotal = f["STotal"];
            SAtally = f["STally"];
            SMatrix = f["SMatrix"];
            #p1Max = f["p3Max"];
            #t1MinMax = f["t3MinMax"];
            close(f)
        else
            SAtotal = zeros(Float64,p1_num,u1_num,p2_num,u2_num);  
            SAtally = zeros(UInt32,p1_num,u1_num,p2_num,u2_num);
            SMatrix = zeros(Float64,p1_num,u1_num,p2_num,u2_num);
            #pMax = zeros(Float64,p1_num,u1_num,p2_num,u2_num);
            #tMinMax = zeros(Float64,2,p1_num,u1_num,p2_num,u2_num);
            #fill!(@view(tMinMax[1,:,:,:,:]),1e0);
        end =#

    # ================================= #

    # ===== Run MonteCarlo Integration ==== #

        prog = Progress(numTiterPerThread)

        if nThreads == 1
            scale = 1.0
        else
            #scale = [range(0.0,1.0,length=10);]
            scale = [0.0;]
        end

for (ii,scale_val) in enumerate(scale)
numT = round(Int,numTiterPerThread/length(scale))

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
        workers  = [EmissionMonteCarloAxi_MultiThread!(GainTotal2,GainTally2,GainTotal3,GainTally3,LossTotal1,LossTally1,ArrayOfLocks,Parameters,numT,numSiterPerThread,nThreads,scale_val,prog,thread) for thread in 1:nThreads]

        wait.(workers) # Allow all workers to finish

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
        u3val = bounds(u_low,u_up,u3_num,u3_grid).*pi
        h3val = bounds(h_low,h_up,h3_num,h3_grid).*pi

        # Apply Momentum space volume elements
        MomentumSpaceFactorsEmission!(LossMatrix1,GainMatrix2,GainMatrix3,p3val,u3val,h3val)

        println("Weighting average of New and Old Sampling Arrays")

        # old arrays are modified in this process 
        WeightedAverageGainEmission!(GainMatrix2,OldGainMatrix2,GainTally2,OldGainTally2,GainMatrix3,OldGainMatrix3,GainTally3,OldGainTally3)
        WeightedAverageLossEmission!(LossMatrix1,OldLossMatrix1,LossTally1,OldLossTally1)

        # preallocate
        #SMatrixOld = copy(SMatrix);
       # fill!(SMatrix,0e0);
        
        # divide element wise by tallys
        #@. SMatrix = SAtotal / SAtally
        #replace!(SMatrix,NaN=>0e0); # remove NaN caused by /0e0

        # Momentum space volume elements and symmetries
        #PhaseSpaceFactorsSync1!(SMatrix,p1val,u1val,p2val,u2val)      # applies phase space factors for symmetries                  
           # initial states are symmetric -> apply symmetry of interaction to improve MC values
        #PhaseSpaceFactorsSync2!(SMatrix,p1val,u1val)            # corrects phase space factors for application in kinetic models
                                            
        # correction to better conserve particle number and account for statistical noise of MC method
        #SCorrection2!(SMatrix,TMatrix) 

        # output a measure of convergence, i.e. new-old/old
        #SConverge = (SMatrix .- SMatrixOld)./SMatrixOld

    # ===================================== #

end # scale loop
finish!(prog)

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