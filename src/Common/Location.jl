"""
    location(low_bound,up_bound,num,val,spacing)

Returns the index of the bin in which 'val' is contained based the grid bounds of that variable with 'num' bins.

Implemented grid spacing types are:
    - uniform spacing: `spacing = "u"`
    - log10 spacing: `spacing = "l"`
        - up and low bounds should be supplied as log10 values
    - binary (1/2^n) spacing: `spacing = "b"` 
        - binary spacing is used for u=cos(theta) grids and therefore bounds should always be [-1 1] and num must be odd!

# Examples
```julia-repl
julia> location(0e0,10e0,9,2e0,"u")
2
```
"""
function location(low_bound::Float64,up_bound::Float64,num::Int64,val::Float64,spacing::String)
    # function for generating position in array
    if spacing == "u" # uniform spacing
        return val != up_bound ? floor(Int64,num*(val-low_bound)/(up_bound-low_bound)+1) : num 
    elseif spacing == "l" # log spacing
        logval = log10(val)
        loc = logval != up_bound ? floor(Int64,num*(logval-low_bound)/(up_bound-low_bound)+1) : num
        return loc 
    elseif spacing == "b" # binary (2^n) fractional spacing
        logval = log(1/2,1-abs(val))
        num_half = Int64((num-1)/2)
        loc = logval < num_half ? floor(Int64,logval/up_bound) : num_half
        return sign(val) == -1 ? num_half+1-loc : num_half+1+loc
    else
        error("Spacing type not recognized")
    end
end

function location(low_bound::Float64,up_bound::Float64,num::Int64,val::Float64,::UniformGrid)
    # grid location for uniform grid
    # if val is on grid boundary then it is assigned to the next bin 
    # if val == up_bound then it is assigned to the last bin
    return val != up_bound ? floor(Int64,num*(val-low_bound)/(up_bound-low_bound)+1) : num
end

#=function location(low_bound::Float64,up_bound::Float64,num::Int64,val::Float64,::LogTenGrid)
    # grid location for log10 grid
    logval = log10(val)
    loc = logval != low_bound ? ceil(Int64,num*(logval-low_bound)/(up_bound-low_bound)) : Int64(1) 
    return 1 <= loc <= num ? loc : loc>num ? num+1 : 1 # assigns 1 for under, num+1 for over and loc for in range
end=#

function location(low_bound::Float64,up_bound::Float64,num::Int64,val::Float64,::LogTenGrid)
    # grid location for log10 grid with no underflow or overflow
    # if val is on grid boundary then it is assigned to the next bin 
    # if val == up_bound then it is assigned to the last bin
    logval::Float64 = log10(val)
    loc::Int64 = logval != up_bound ? floor(Int64,num*(logval-low_bound)/(up_bound-low_bound)+1) : num
    return loc 
end

function location(low_bound::Float64,up_bound::Float64,num::Int64,val::Float64,::BinaryGrid)
    # grid location for binary grid
    logval = log(1/2,1-abs(val))
    num_half = Int64((num-1)/2)
    loc = logval < num_half ? floor(Int64,logval/up_bound) : num_half
    return sign(val) == -1 ? Int64(num_half+1-loc) : Int64(num_half+1+loc)
end

function Grid_String_to_Type(grid_string)
    if grid_string == "u"
        return UniformGrid()
    elseif grid_string == "l"
        return LogTenGrid()
    elseif grid_string == "b"
        return BinaryGrid()
    else
        error("Spacing type not recognized")
    end
end

"""
    location_t(numt,val)

Returns the index of the bin in which the costheta 'val' is contatined based on the 'numt' of bins. Bounds [tl tu] are defined as CONST in Init.jl

# Examples
```julia-repl
julia> location_t(8,0.5e0)
6
```
"""
function location_t(numt::Int64,val::Float64)
    # function for generating position in array. Bins MUST be uniform
    return val != u_low ? ceil(Int64,Float64(numt)*(val-u_low)/(u_up-u_low)) : Int64(1) 
end

"""
    location_p(u,l,num,val)

Returns the index of the momentum bin in which 'val' is contatined based on the 'num' of bins and their 'u' upper and 'l' lower bound including overflow and underflow possibilities. Overflow are assigned to num+1 while underflow are assigned to lowest bin i.e. 1.

# Examples
```julia-repl
julia> location_p(10e0,1e0,9,2e0)
2
julia> location_p(10e0,1e0,9,11e0) # overflow
10
julia> location_p(10e0,1e0,9,0.5e0) # underflow
1
```
"""
function location_p(u::Float64,l::Float64,num::Int64,val::Float64)
    # function for generating poisition in array. Bins MUST be uniform
    logp = log10(val)
    loc = logp != l ? ceil(Int64,Float64(num)*(logp-l)/(u-l)) : Int64(1) 
    return 1 <= loc <= num ? loc : loc>num ? num+1 : 1 # assignes 1 for under, num+1 for over and loc for in range
end

"""
    vectorLocation(pu,pl,nump,numt,vector)

Returns a tuple of bin location for (log10momentum,cos(theta)) based on an input 'vector' and bounds 'u,l' of their domains and the 'num' of uniformly spaced bins.
costheta bounds [tl tu] are defined as CONST in Init.jl

# Examples
```julia-repl
julia> vectorLocation(4e0,-5e0,9,8,[1e0,0.5e0,1.5e0])
(5,6)
"""
function vectorLocation(pu::Float64,pl::Float64,nump::Int64,numt::Int64,vector::Vector{Float64})

    logp = log10(vector[1])
    ctheta = vector[2]

    logploc = (logp != pl ? ceil(Int64,Float64(nump)*(logp-pl)/(pu-pl)) : Int64(1))
    cthetaloc = (ctheta != tl ? ceil(Int64,Float64(numt)*(ctheta-tl)/(tu-tl)) : Int64(1))
    
    return (logploc,cthetaloc)
    
end


