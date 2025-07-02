function UserBinaryParameters()

    #name1::String = eval(Symbol("name1"))

    name1::String = getfield(Main,Symbol("name1"))
    name2::String = getfield(Main,Symbol("name2"))
    name3::String = getfield(Main,Symbol("name3"))
    name4::String = getfield(Main,Symbol("name4"))

    p1_low::Float64 = Symbol("p_low_"*name1)
    p1_up::Float64 = getfield(Main,Symbol("p_up_"*name1))
    p1_grid::String = getfield(Main,Symbol("p_grid_"*name1))
    p1_num::Int64 = getfield(Main,Symbol("p_num_"*name1))

    u1_grid::String = getfield(Main,Symbol("u_grid_"*name1))
    u1_num::Int64 = getfield(Main,Symbol("u_num_"*name1))
    
    h1_grid::String = getfield(Main,Symbol("h_grid_"*name1))
    h1_num::Int64 = getfield(Main,Symbol("h_num_"*name1))

    p2_low::Float64 = getfield(Main,Symbol("p_low_"*name2))
    p2_up::Float64 = getfield(Main,Symbol("p_up_"*name2))
    p2_grid::String = getfield(Main,Symbol("p_grid_"*name2))
    p2_num::Int64 = getfield(Main,Symbol("p_num_"*name2))

    u2_grid::String = getfield(Main,Symbol("u_grid_"*name2))
    u2_num::Int64 = getfield(Main,Symbol("u_num_"*name2))

    h2_grid::String = getfield(Main,Symbol("h_grid_"*name2))
    h2_num::Int64 = getfield(Main,Symbol("h_num_"*name2))

    p3_low::Float64 = getfield(Main,Symbol("p_low_"*name3))
    p3_up::Float64 = getfield(Main,Symbol("p_up_"*name3))
    p3_grid::String = getfield(Main,Symbol("p_grid_"*name3))
    p3_num::Int64 = getfield(Main,Symbol("p_num_"*name3))

    u3_grid::String = getfield(Main,Symbol("u_grid_"*name3))
    u3_num::Int64 = getfield(Main,Symbol("u_num_"*name3))

    h3_grid::String = getfield(Main,Symbol("h_grid_"*name3))
    h3_num::Int64 = getfield(Main,Symbol("h_num_"*name3))

    p4_low::Float64 = getfield(Main,Symbol("p_low_"*name4))
    p4_up::Float64 = getfield(Main,Symbol("p_up_"*name4))
    p4_grid::String = getfield(Main,Symbol("p_grid_"*name4))
    p4_num::Int64 = getfield(Main,Symbol("p_num_"*name4))

    u4_grid::String = getfield(Main,Symbol("u_grid_"*name4))
    u4_num::Int64 = getfield(Main,Symbol("u_num_"*name4))

    h4_grid::String = getfield(Main,Symbol("h_grid_"*name4))
    h4_num::Int64 = getfield(Main,Symbol("h_num_"*name4))

    m1::Float64 = getfield(DiplodocusCollisions,Symbol("mu"*name1))
    m2::Float64 = getfield(DiplodocusCollisions,Symbol("mu"*name2))
    m3::Float64 = getfield(DiplodocusCollisions,Symbol("mu"*name3))
    m4::Float64 = getfield(DiplodocusCollisions,Symbol("mu"*name4))

    Parameters = (name1,name2,name3,name4,m1,m2,m3,m4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num,h4_grid,h4_num)

    numLoss::Int64 = getfield(Main,Symbol("numLoss"))   
    numGain::Int64 = getfield(Main,Symbol("numGain"))  
    numThreads::Int64 = getfield(Main,Symbol("numThreads"))

    fileLocation::String = getfield(Main,Symbol("fileLocation"))
    fileName::String = FileName(Parameters)

    Setup::Tuple = (Parameters, numLoss, numGain, numThreads, fileLocation, fileName)

    return Setup
end

"""
    FileName(Parameters)

Generates a file name based on user provided parameters
"""
function FileName(Parameters::Tuple{String,String,String,String,Float64,Float64,Float64,Float64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64})
    
    (name1,name2,name3,name4,mu1,mu2,mu3,mu4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num,h4_grid,h4_num) = Parameters
    
    fileName = name1*name2*name3*name4
    fileName *= "#"*string(p1_low)*"-"*string(p1_up)*p1_grid*string(p1_num)
    fileName *= "#"*u1_grid*string(u1_num)
    fileName *= "#"*h1_grid*string(h1_num)

    fileName *= "#"*string(p2_low)*"-"*string(p2_up)*p2_grid*string(p2_num)
    fileName *= "#"*u2_grid*string(u2_num)
    fileName *= "#"*h2_grid*string(h2_num)

    fileName *= "#"*string(p3_low)*"-"*string(p3_up)*p3_grid*string(p3_num)
    fileName *= "#"*u3_grid*string(u3_num)
    fileName *= "#"*h3_grid*string(h3_num)

    fileName *= "#"*string(p4_low)*"-"*string(p4_up)*p4_grid*string(p4_num)
    fileName *= "#"*u4_grid*string(u4_num)
    fileName *= "#"*h4_grid*string(h4_num)
    
    fileName *= ".jld2";

    return fileName
end