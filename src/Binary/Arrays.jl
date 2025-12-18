"""
    MonteCarloArraysBinary(Parameters)

Generates arrays for Monte Carlo sampling for binary (12->34) interactions.
"""
function MonteCarloArraysBinary(Parameters::Tuple{String,String,String,String,Float64,Float64,Float64,Float64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64})

    (name1,name2,name3,name4,mu1,mu2,mu3,mu4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num,h4_grid,h4_num) = Parameters

    GainTotal3::Array{Float64,9} = zeros(Float64,p3_num+2,u3_num,h3_num,p1_num,u1_num,h1_num,p2_num,u2_num,h2_num); 
    LossTotal::Array{Float64,6} = zeros(Float64,p1_num,u1_num,h1_num,p2_num,u2_num,h2_num);
    # Gain arrays have first dimension elements [k_underflow,k1,k2,k3,...,kn,k_overflow,N]
    GainTally3::Array{UInt32,9} = zeros(UInt32,(p3_num+3),u3_num,h3_num,p1_num,u1_num,h1_num,p2_num,u2_num,h2_num);
    LossTally::Array{UInt32,6} = zeros(UInt32,p1_num,u1_num,h1_num,p2_num,u2_num,h2_num);

    GainMatrix3::Array{Float64,9} = zeros(Float64,p3_num+2,u3_num,h3_num,p1_num,u1_num,h1_num,p2_num,u2_num,h2_num);
    LossMatrix1::Array{Float64,6} = zeros(Float64,p1_num,u1_num,h1_num,p2_num,u2_num,h2_num);
    LossMatrix2::Array{Float64,6} = zeros(Float64,p2_num,u2_num,h2_num,p1_num,u1_num,h1_num);

    GainTotal4::Array{Float64,9} = (mu3 == mu4) ? zeros(Float64,0,0,0,0,0,0,0,0,0) : zeros(Float64,p4_num+2,u4_num,h4_num,p1_num,u1_num,h1_num,p2_num,u2_num,h2_num); 
    GainTally4::Array{UInt32,9} = (mu3 == mu4) ? zeros(UInt32,0,0,0,0,0,0,0,0,0) : zeros(UInt32,(p4_num+3),u4_num,h4_num,p1_num,u1_num,h1_num,p2_num,u2_num,h2_num);
    GainMatrix4::Array{Float64,9} = (mu3 == mu4) ? zeros(Float64,0,0,0,0,0,0,0,0,0) : zeros(Float64,p4_num+2,u4_num,h4_num,p1_num,u1_num,h1_num,p2_num,u2_num,h2_num);

    return (GainTotal3,GainTotal4,LossTotal,GainTally3,GainTally4,LossTally,GainMatrix3,GainMatrix4,LossMatrix1,LossMatrix2)

end


"""
    OldMonteCarloArraysBinary(Parameters)

Load/Generates arrays for Monte Carlo sampling for binary (12->34) interactions for if there is/is not a previously saved sample.
"""
function OldMonteCarloArraysBinary(Parameters::Tuple{String,String,String,String,Float64,Float64,Float64,Float64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64},filePath::String)

    (name1,name2,name3,name4,mu1,mu2,mu3,mu4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num,h4_grid,h4_num) = Parameters
    fileExist = isfile(filePath)

    if fileExist
        f = jldopen(filePath,"r+");
        OldGainWeights3 = f["GainWeights3"];
        OldGainMatrix3 = f["GainMatrix3"];
        OldGainWeights4 = f["GainWeights4"];
        OldGainMatrix4 = f["GainMatrix4"];
        OldLossTally = f["LossTally"];
        OldLossMatrix1 = f["LossMatrix1"];
        OldLossMatrix2 = f["LossMatrix2"];
        close(f)
    else
        # Gain arrays have first dimension elements [w_underflow,w1,w2,w3,...,wn,w_overflow]
        OldGainWeights3::Array{Float64,9} = zeros(UInt32,p3_num+2,u3_num,h3_num,p1_num,u1_num,h1_num,p2_num,u2_num,h2_num);
        OldGainMatrix3::Array{Float64,9} = zeros(Float64,p3_num+2,u3_num,h3_num,p1_num,u1_num,h1_num,p2_num,u2_num,h2_num);
        OldGainWeights4::Array{Float64,9} = zeros(UInt32,p4_num+2,u4_num,h4_num,p1_num,u1_num,h1_num,p2_num,u2_num,h2_num);
        OldGainMatrix4::Array{Float64,9} = zeros(Float64,p4_num+2,u4_num,h4_num,p1_num,u1_num,h1_num,p2_num,u2_num,h2_num);

        OldLossTally::Array{UInt32,6} = zeros(UInt32,p1_num,u1_num,h1_num,p2_num,u2_num,h2_num);
        OldLossMatrix1::Array{Float64,6} = zeros(Float64,p1_num,u1_num,h1_num,p2_num,u2_num,h2_num);
        OldLossMatrix2::Array{Float64,6} = zeros(Float64,p2_num,u2_num,h2_num,p1_num,u1_num,h1_num);
    end

    return (OldGainWeights3,OldGainWeights4,OldLossTally,OldGainMatrix3,OldGainMatrix4,OldLossMatrix1,OldLossMatrix2)

end