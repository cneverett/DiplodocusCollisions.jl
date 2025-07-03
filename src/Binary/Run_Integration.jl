#=
This script is provided as an example of how to run the BoltzmannCollisionIntegral.jl module.
=#

# First select the binary interaction (12->34) you want to evaluate the collision integral for by setting the names of particles 1 through 4. The pairs of names (12) and (34) should be in alphabetical order.

    name1::String = "Sph";
    name2::String = "Sph";
    name3::String = "Sph";
    name4::String = "Sph";

# Define the Momentum space discretisation for each particle species named. This includes the upper and lower momentum bounds, grid type "l", "u", or "b", and the number of bins. Must be of the format `p_low_name`, p_low_name`, `p_grid_name`, `p_num_name`, `u_grid_name` and `u_num_name` where `name` is the abbreviated three letter name of the particle species. 

    p_low_Sph = -5e0
    p_up_Sph = 4e0
    p_grid_Sph = "l"
    p_num_Sph = 72

    u_grid_Sph = "u"
    u_num_Sph = 8

    phi_grid_Sph = "u"
    phi_num_Sph = 1

# ==== DO NOT EDIT THIS SECTION ======= #
# ===================================== #

    # This section takes all the information above, and generates a tuple of the parameters to be passed to the evaluation functions.

    p1_low::Float64 = getfield(Main,Symbol("p_low_"*name1))
    p1_up::Float64 = getfield(Main,Symbol("p_up_"*name1))
    p1_grid::String = getfield(Main,Symbol("p_grid_"*name1))
    p1_num::Int64 = getfield(Main,Symbol("p_num_"*name1))

    u1_grid::String = getfield(Main,Symbol("u_grid_"*name1))
    u1_num::Int64 = getfield(Main,Symbol("u_num_"*name1))
    
    h1_grid::String = getfield(Main,Symbol("phi_grid_"*name1))
    h1_num::Int64 = getfield(Main,Symbol("phi_num_"*name1))

    p2_low::Float64 = getfield(Main,Symbol("p_low_"*name2))
    p2_up::Float64 = getfield(Main,Symbol("p_up_"*name2))
    p2_grid::String = getfield(Main,Symbol("p_grid_"*name2))
    p2_num::Int64 = getfield(Main,Symbol("p_num_"*name2))

    u2_grid::String = getfield(Main,Symbol("u_grid_"*name2))
    u2_num::Int64 = getfield(Main,Symbol("u_num_"*name2))

    h2_grid::String = getfield(Main,Symbol("phi_grid_"*name2))
    h2_num::Int64 = getfield(Main,Symbol("phi_num_"*name2))

    p3_low::Float64 = getfield(Main,Symbol("p_low_"*name3))
    p3_up::Float64 = getfield(Main,Symbol("p_up_"*name3))
    p3_grid::String = getfield(Main,Symbol("p_grid_"*name3))
    p3_num::Int64 = getfield(Main,Symbol("p_num_"*name3))

    u3_grid::String = getfield(Main,Symbol("u_grid_"*name3))
    u3_num::Int64 = getfield(Main,Symbol("u_num_"*name3))

    h3_grid::String = getfield(Main,Symbol("phi_grid_"*name3))
    h3_num::Int64 = getfield(Main,Symbol("phi_num_"*name3))

    p4_low::Float64 = getfield(Main,Symbol("p_low_"*name4))
    p4_up::Float64 = getfield(Main,Symbol("p_up_"*name4))
    p4_grid::String = getfield(Main,Symbol("p_grid_"*name4))
    p4_num::Int64 = getfield(Main,Symbol("p_num_"*name4))

    u4_grid::String = getfield(Main,Symbol("u_grid_"*name4))
    u4_num::Int64 = getfield(Main,Symbol("u_num_"*name4))

    h4_grid::String = getfield(Main,Symbol("phi_grid_"*name4))
    h4_num::Int64 = getfield(Main,Symbol("phi_num_"*name4))

    mu1::Float64 = getfield(BoltzmannCollisionIntegral,Symbol("mu"*name1))
    mu2::Float64 = getfield(BoltzmannCollisionIntegral,Symbol("mu"*name2))
    mu3::Float64 = getfield(BoltzmannCollisionIntegral,Symbol("mu"*name3))
    mu4::Float64 = getfield(BoltzmannCollisionIntegral,Symbol("mu"*name4))

    Parameters = BinaryParameters(name1,name2,name3,name4,mu1,mu2,mu3,mu4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num,h4_grid,h4_num)

# ===================================== #
# ===================================== #

# Define the number of Monte-Carlo samples to take. numTiter gives the number of random sets of (p1,p2) vectors that are sampled, while numSiter given the number of random (p3) vectors that are sampled for each set of (p1,p2) vectors. For good convergence, these values should be much larger than the number of bins in the momentum space discretisation. 
# As a rule of thumb, each single iteration takes approx 250ns to complete on a single core of a modern CPU

    # For Serial Evaluation
    numTiter::Int64 = p1_num*u1_num*h1_num*p2_num*u2_num*h2_num*1;    # number of T matrix iterations i.e. random p1 p2 points
    numSiter::Int64 = p3_num*u3_num*h3_num*1;        # number of S matrix iteration per T matrix iteration i.e. random p3 directions per p1 p2 point

    # For MultiThread Evaluation
    numTiterPerThread::Int64 = p1_num*u1_num*h1_num*p2_num*u2_num*h2_num*1;;    # number of T matrix iterations i.e. random p1 p2 points. Should be > nump1*numt1*nump2*numt2 to ensure good sampling.
    numSiterPerThread::Int64 = p3_num*u3_num*h3_num*1;        # number of S matrix iteration per T matrix iteration i.e. random p3 directions per p1 p2 point. Should be > nump3*numt3 to ensure good sampling.

# Define the number of Threads (more specifically, workers) used for MultiThreaded evaluation. (see https://docs.julialang.org/en/v1/manual/multi-threading/ for how to set up multi-threading in Julia)

    nThreads::Int64 = 10;

# Define the location where the output file (.jld2 format) is to be stored and its name.
# By default the file name contains all the information of discretisation in form p1_low-p1_up p1_grid p1_num#...#p4_low-p4_up p4_grid nump4#u1_grid u1_num#...#u4_grid u4_num, and it to be stored in a folder called "Data" in the current working directory.

    fileLocation::String = pwd()*"\\Data";
    #fileName::String = name1*name2*name3*name4*"#"*string(p1_low)*"-"*string(p1_up)*p1_grid*string(p1_num)*"#"*string(p2_low)*"-"*string(p2_up)*p2_grid*string(p2_num)*"#"*string(p3_low)*"-"*string(p3_up)*p3_grid*string(p3_num)*"#"*string(p4_low)*"-"*string(p4_up)*p4_grid*string(p4_num)*"#"*u1_grid*string(u1_num)*"#"*u2_grid*string(u2_num)*"#"*u3_grid*string(u3_num)*"#"*u4_grid*string(u4_num)*".jld2";
    MinMax = false;

# Now run the evaluation functions.

    # ===== DO NOT EDIT THIS SECTION ===== #
    
        userInputSerial = BinaryUserInput(Parameters,numTiter,numSiter,fileLocation;MinMax=MinMax)
        userInputMultiThread = BinaryUserInput(Parameters,numTiterPerThread,numSiterPerThread,fileLocation;nThreads=nThreads,MinMax=MinMax)

    # ===================================== #

    # For Serial Evaluation
    SpectraEvaluateSerial(userInputSerial)

    # For MultiThread Evaluation
    #SpectraEvaluateMultiThread(userInputMultiThread)

# Once evaluation is complete all the data stored in the output file can be loaded into the user workspace using the following function:

    #(Stot,Ttot,Stal,Ttal,SMatrix,TMatrix,p3Max,t3MinMax,SConv,TConv) = fload_All(fileLocation,fileName);


