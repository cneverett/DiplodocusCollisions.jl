#= Functions for the S and T integration functions =#

"""
    CosThetaValue(p1v,p2v)

Returns the cosine of the angle between two momentum vectors `p1v` and `p2v` of format [p,cos(theta),phi/pi,theta/pi]. To aid in floating point precision, 1 is subtracted from the returned cosine value. i.e. CosTheta12m1 = cos(Theta12) - 1.0. This is done as it is often the case that two momentum states are very close together.
"""
function ThetaValue(p1v::Vector{Float64},p2v::Vector{Float64})

    t1::Float64 = p1v[4]
    h1::Float64 = p1v[3]
    t2::Float64 = p2v[4]
    h2::Float64 = p2v[3]
    
    if (td::Float64=pi*abs(t1-t2)/2) < 1e-6 && (hd::Float64=pi*abs(h1-h2)/2) < 1e-6
        ta::Float64 = (t1 + t2)/2
        CosTheta12m1 = 2*td^2*(hd^2-1)-2*hd^2*sinpi(ta)^2
    else
        CosTheta12m1 = cospi(t1)*cospi(t2) + sinpi(t1)*sinpi(t2)*cospi(h1-h2) - 1e0
    end

    return CosTheta12m1

end

"""
    TValue(p1v,p2v,sigma,mu1,mu2)

returns `Tval` with its Tval from MC integration based on initial momentum states `p1v` and `p2v` and cross section `sigma` based on particle selection.
```math
T_\\text{val} = \\frac{1}{p^0_1p^0_2}\\sigma(s)F_12(s)
```
If initial state fails `sCheck`, i.e. cannot generate a physical output state, Tval is set to 0e0. 
Assumes f(x,p,u,ϕ)=f(x,vec{p})/p^2=constant over bin
"""
function TValue(p1v::Vector{Float64},p2v::Vector{Float64},sigma::Function,m1::Float64,m2::Float64,m3::Float64,m4::Float64)

    # pre-defining terms for efficiency 
    p1::Float64 = p1v[1]
    p2::Float64 = p2v[1]

    st1::Float64,ct1::Float64 = sincospi(p1v[4])
    st2::Float64,ct2::Float64 = sincospi(p2v[4])

    #ct1::Float64 = p1v[2] #cospi(p1v[2])
    #ct2::Float64 = p2v[2] #cospi(p2v[2]) 

    #st1::Float64 = sqrt(1e0-p1v[2]^2) #sinpi(p1v[2])
    #st2::Float64 = sqrt(1e0-p2v[2]^2) #sinpi(p2v[2])

    ch1h2::Float64 = cospi(p1v[3]-p2v[3])

    ctheta12::Float64 = ct1*ct2+ch1h2*st1*st2

    m12::Float64 = m1^2
    m22::Float64 = m2^2

    Tval::Float64 = 0e0
    
    #= This method leads to errors when p/m < 10^-6 due to Floating point precision, causing s<(m+m)^2 and InvFlux to return a complex number error 
    E1 = sqrt(p1^2+m12)
    E2 = sqrt(p2^2+m22) 
    s = m12+m22 + 2*E1*E2 - 2*p1*p2*(ct1*ct2+ch1h2*st1*st2)
    =#

    #= will attempt to use the algebraic relation sqrt(1+x)-1 = x/(sqrt(1+x)+1) to split E up into a large and small part i.e. sqrt(p^2+m^2) = m + (p^2)/(sqrt(m^2+p^2)+m) = E 
    then label Es = (p^2)/(sqrt(m^2+p^2)+m) and rewrite 
    s = (m3+m4)^2 + 2(m3*Es4+m4*Es3+Es3*Es4 - p3p4(...)) =#

    Es1::Float64 = (p1^2)/(sqrt(m12+p1^2)+m1)
    Es1s::Float64 = Es1/p1
    Es2::Float64 = (p2^2)/(sqrt(m22+p2^2)+m2)
    Es2s::Float64 = Es2/p2

    sBig::Float64 = (m1+m2)^2
    sSmol::Float64 = 2*p1*p2*(-ctheta12 +Es1s*Es2s +m1*Es2s/p1 +m2*Es1s/p2)

    if sCheck(sSmol,sBig,m1,m2,m3,m4) # check if s value is valid for interaction

        E1::Float64 = Es1 + m1
        E2::Float64 = Es2 + m2
        
        Tval = (1/E1)*(1/E2)*(InvariantFluxSmall(sSmol,m1,m2))*sigma(sSmol,sBig)
        if Tval==Inf || Tval < 0e0
            println("")
            println("p1v = $p1v")
            println("p2v = $p2v")
            println("Tval = $Tval")
            println("sSmol = $sSmol")
            println("sBig = $sBig")
            error("ST1 Inf or -ve#")
        end
    else # if not valid set T value to zero
        Tval = 0e0
    end

    return Tval, sBig, sSmol

end


"""
    SValue3(p3v,p1v,p2v,dsigmadt,mu1,mu2,mu3)

Returns `Sval` from MC integration based on initial momentum states `p1v` and `p2v` and final state `p3v` and differential cross section `dsigmadt` based on particle selection 12->34. 
```math
S_\\text{val}=\\frac{\\mathrm{d}\\sigma_{12|34}}{\\mathrm{d}t}\\frac{\\mathcal{F}_{12}^2}{\\pi\\left|p^0_3(p_1\\cos\\Theta_{31}+p_2\\cos\\Theta_{32})-p_3(p^0_1+p^0_2)\\right|}\frac{p_3^2}{p^0_1p^0_2}.
``` 
Assumes f(x,p,u,ϕ)=f(x,vec{p})/p^2=constant over bin
"""
function SValue3(p3v::Vector{Float64},p1v::Vector{Float64},p2v::Vector{Float64},sBig::Float64,sSmol::Float64,dsigmadt::Function,m1::Float64,m2::Float64,m3::Float64,m4::Float64,prob::Float64)

    # pre-defining terms for efficiency 
    p1::Float64 = p1v[1]
    p2::Float64 = p2v[1]
    st1::Float64,ct1::Float64 = sincospi(p1v[4])
    st2::Float64,ct2::Float64 = sincospi(p2v[4])

    #ct1::Float64 = p1v[2] #cospi(p1v[2])
    #ct2::Float64 = p2v[2] #cospi(p2v[2]) 

    #st1::Float64 = sqrt(1e0-p1v[2]^2) #sinpi(p1v[2])
    #st2::Float64 = sqrt(1e0-p2v[2]^2) #sinpi(p2v[2])

    ch1h2::Float64 = cospi(p1v[3]-p2v[3])

    p3::Float64 = p3v[1]
    (st3::Float64,ct3::Float64) = sincospi(p3v[4])
    ch3h1::Float64 = cospi(p3v[3]-p1v[3])
    ch3h2::Float64 = cospi(p3v[3]-p2v[3])

    m32 = m3^2
    m42 = m4^2
    m12 = m1^2
    m22 = m2^2
    
    Es1::Float64 = m1 != 0e0 ? (p1^2)/(sqrt(m12+p1^2)+m1) : p1
    Es1s::Float64 = Es1/p1
    E1::Float64 = Es1 + m1
 
    Es2::Float64 = m2 != 0e0 ? (p2^2)/(sqrt(m22+p2^2)+m2) : p2
    Es2s::Float64 = Es2/p2
    E2::Float64 = Es2 + m2

    Es3::Float64 = m3 != 0e0 ? (p3^2)/(sqrt(m32+p3^2)+m3) : p3
    Es3s::Float64 = Es3/p3
    E3::Float64 = Es3 + m3

    ctheta12::Float64 = ct1*ct2+ch1h2*st1*st2
    ctheta13::Float64 = ct3*ct1+ch3h1*st3*st1
    ctheta23::Float64 = ct3*ct2+ch3h2*st3*st2
    
    deltacorrect::Float64 = Es1*p3 - Es3*p1*ctheta13
    deltacorrect += m1*p3 - m3*p1*ctheta13
    deltacorrect += Es2*p3 - Es3*p2*ctheta23    
    deltacorrect += m2*p3 - m3*p2*ctheta23

    #sBig::Float64 = (m1+m2)^2
    #sSmol::Float64 = 2*(m1*Es2 + m2*Es1 + Es1*Es2 - p1*p2*(ct1*ct2+ch1h2*st1*st2))
    #sSmol::Float64 = 2*p1*p2*(-ctheta12 + Es1s*Es2s + m1*Es2s/p1 + m2*Es1s/p2)

    # t = tBig + tSmol
    tBig::Float64 = (m3-m1)^2
    #tSmol::Float64 = -2*(m1*Es3 + m3*Es1 + Es3*Es1 - p3*p1*(ct3*ct1+ch3h1*st3*st1))
    #tSmol::Float64 = -2*m1*m3 - 2*E1*E3 + 2*p1*p3*ctheta13
    tSmol::Float64 = 2*p3*p1*(ctheta13 - Es3s*Es1s - m1*Es3s/p1 - m3*Es1s/p3)
    # u = uBig + uSmol
    uBig::Float64 = (m2-m3)^2
    #uSmol::Float64 = m12+m22+m32+m42 - sBig - tBig - uBig - sSmol - tSmol  # this leads to Float64 issues better to calculate directly
    #uSmol::Float64 = -2*(m3*Es2 + m2*Es3 + Es2*Es3 - p2*p3*(ct2*ct3+ch3h2*st2*st3))
    uSmol::Float64 = 2*p2*p3*(ctheta23 - Es2s*Es3s - m3*Es2s/p3 - m2*Es3s/p2)

    # Sspe anisotropic emission spectrum (to be integrated over d^2p1d^3p3d^3p4). See obsidian note on discrete anisotropic kinetic equation
    val::Float64 = (1/E1)*(1/E2)*(InvariantFlux2Small(sSmol,m1,m2))/pi

    #deltacorrect::Float64 = (Es1*p3 - Es3*p1*(ct3*ct1+ch3h1*st3*st1) - m3*p1*(ct3*ct1+ch3h1*st3*st1)) + (Es2*p3 - Es3*p2*(ct3*ct2+ch3h2*st3*st2) - m3*p2*(ct3*ct2+ch3h2*st3*st2)) + m1*p3 + m2*p3
    # more Float accurate for when p1 and p2 have large order of magnitude difference as sum uses pairwise summation to reduce round of errors
    #sumTerms .= (Es1, -Es3prime*p1*(ct3*ct1+ch3h1*st3*st1), Es2, -Es3prime*p2*(ct3*ct2+ch3h2*st3*st2), m1, -(m3/p3)*p1*(ct3*ct1+ch3h1*st3*st1), m2, -(m3/p3)*p2*(ct3*ct2+ch3h2*st3*st2))
    #deltacorrect = p3*sum_oro(sumTerms)

    #=if (stuCheck(sSmol,sBig,tSmol,tBig,uSmol,uBig,m1,m2,m3,m4) == false || tCheck(tSmol,tBig,m1,m2,m3,m4) == false || uCheck(uSmol,uBig,m1,m2,m3,m4) == false)
        println("p1v = $p1v")
        println("p2v = $p2v")
        println("p3v = $p3v")
        error("stu check")
    end=#

    Sval = dsigmadt(sSmol,sBig,tSmol,tBig,uSmol,uBig)*val*(p3^2/(deltacorrect*sign(deltacorrect)))

    #println("S3:",Sval,p3v,p1v,p2v)
    #println(sSmol)
    #println(tSmol)
    #println(uSmol)
    #println(sSmol+uSmol+tSmol)
    #println(string(val)*"#",string(p3^2)*"#",string(1/deltacorrect)*"#",string(dsigmadt(sSmol,sBig,tSmol,tBig,uSmol,uBig))*"#")

    if (Sval==Inf || Sval == -Inf || Sval < 0e0)
        println("")
        println("Sval = $Sval")
        println("deltacorrect = $deltacorrect")
        println("dsdt = $(dsigmadt(sSmol,sBig,tSmol,tBig,uSmol,uBig))")
        println("sSmol = $sSmol")
        println("tSmol = $tSmol")
        println("uSmol = $uSmol")
        println("p1v = $p1v")
        println("p2v = $p2v")
        println("p3v = $p3v")
        error("Sval Inf")  
    end

    #=     if (Sval/prob) >= 1e13
        println("")
        println("sSmol = $sSmol")
        println("tSmol = $tSmol")
        println("uSmol = $uSmol")
        println("check = $(sSmol+tSmol+uSmol+sBig+tBig+uBig)")
    end =#

    return Sval

end


"""
    SValue4(p3v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)

Returns `Sval` from MC integration based on initial momentum states `p1v` and `p2v` and final state `p4v` and differential cross section `dsigmadt` based on particle selection 12->34.  
Assumes f(x,p,μ)=constant over bin
"""
function SValue4(p4v::Vector{Float64},p1v::Vector{Float64},p2v::Vector{Float64},sBig::Float64,sSmol::Float64,dsigmadt::Function,m1::Float64,m2::Float64,m3::Float64,m4::Float64,prob::Float64)

    # pre-defining terms for efficiency 
    p1::Float64 = p1v[1]
    p2::Float64 = p2v[1]
    st1::Float64,ct1::Float64 = sincospi(p1v[4])
    st2::Float64,ct2::Float64 = sincospi(p2v[4])

    #ct1::Float64 = p1v[2] #cospi(p1v[2])
    #ct2::Float64 = p2v[2] #cospi(p2v[2]) 

    #st1::Float64 = sqrt(1e0-p1v[2]^2) #sinpi(p1v[2])
    #st2::Float64 = sqrt(1e0-p2v[2]^2) #sinpi(p2v[2])

    ch1h2::Float64 = cospi(p1v[3]-p2v[3])

    p4::Float64 = p4v[1]
    st4::Float64,ct4::Float64 = sincospi(p4v[4])
    ch4h1::Float64 = cospi(p4v[3]-p1v[3])
    ch4h2::Float64 = cospi(p4v[3]-p2v[3])

    m32 = m3^2
    m42 = m4^2
    m12 = m1^2
    m22 = m2^2

    Es1::Float64 = m1 != 0e0 ? (p1^2)/(sqrt(m12+p1^2)+m1) : p1
    Es1s::Float64 = Es1/p1
    E1::Float64 = Es1 + m1
 
    Es2::Float64 = m2 != 0e0 ? (p2^2)/(sqrt(m22+p2^2)+m2) : p2
    Es2s::Float64 = Es2/p2
    E2::Float64 = Es2 + m2

    Es4::Float64 = m4 != 0e0 ? (p4^2)/(sqrt(m42+p4^2)+m4) : p4
    Es4s::Float64 = Es4/p4
    E4::Float64 = Es4 + m4

    ctheta12::Float64 = ct1*ct2+ch1h2*st1*st2
    ctheta14::Float64 = ct4*ct1+ch4h1*st4*st1
    ctheta24::Float64 = ct4*ct2+ch4h2*st4*st2

    deltacorrect::Float64 = Es1*p4 - Es4*p1*ctheta14
    deltacorrect += m1*p4 - m4*p1*ctheta14
    deltacorrect += Es2*p4 - Es4*p2*ctheta24
    deltacorrect += m2*p4 - m4*p2*ctheta24

    #sBig::Float64 = (m1+m2)^2
    #sSmol::Float64 = 2*(m1*Es2 + m2*Es1 + Es1*Es2 - p1*p2*(ct1*ct2+ch1h2*st1*st2))
    #sSmol::Float64 = 2*p1*p2*(-ctheta12 + Es1s*Es2s + m1*Es2s/p1 + m2*Es1s/p2)

    # u = uBig + uSmol
    uBig::Float64 = (m4-m1)^2
    #uSmol::Float64 = -2*(m1*Es4 + m4*Es1 + Es4*Es1 - p4*p1*(ct4*ct1+ch4h1*st4*st1))
    uSmol::Float64 = 2*p4*p1*(ctheta14 - Es4s*Es1s - m1*Es4s/p1 - m4*Es1s/p4)
    # t = tBig + tSmol
    tBig::Float64 = (m2-m4)^2
    #tSmol::Float64 = m12+m22+m32+m42 - sBig - uBig - tBig - sSmol - uSmol # Leads to Float64 issues, better to calculate directly 
    #tSmol::Float64 = -2*(m4*Es2 + m2*Es4 + Es2*Es4 - p2*p4*(ct2*ct4+ch4h2*st2*st4))
    tSmol::Float64 = 2*p2*p4*(ctheta24 - Es2s*Es4s - m4*Es2s/p4 - m2*Es4s/p2)
    
    # Sspe anisotropic emission spectrum (to be integrated over d^2p1d^3p3d^3p4). See obsidian note on discrete anisotropic kinetic equation
    val::Float64 = (1/E1)*(1/E2)*(InvariantFlux2Small(sSmol,m1,m2))/pi

    # more Float accurate for when p1 and p2 have large order of magnitude difference as sum uses pairwise summation to reduce round of errors
    #sumTerms .= (Es1, -Es3prime*p1*(ct3*ct1+ch3h1*st3*st1), Es2, -Es3prime*p2*(ct3*ct2+ch3h2*st3*st2), m1, -(m3/p3)*p1*(ct3*ct1+ch3h1*st3*st1), m2, -(m3/p3)*p2*(ct3*ct2+ch3h2*st3*st2))
    #deltacorrect = p3*sum_oro(sumTerms)

    #=if (stuCheck(sSmol,sBig,tSmol,tBig,uSmol,uBig,m1,m2,m3,m4) == false || tCheck(tSmol,tBig,m1,m2,m3,m4) == false || uCheck(uSmol,uBig,m1,m2,m3,m4) == false)
        error("stu check")
    end=#

    Sval = dsigmadt(sSmol,sBig,tSmol,tBig,uSmol,uBig)*val*(p4^2/(deltacorrect*sign(deltacorrect)))

    #println("S4:",Sval,p4v,p1v,p2v)
    #println(sSmol)
    #println(tSmol)
    #println(uSmol)
    #println(sSmol+uSmol+tSmol)
    #println(string(val)*"#",string(p4^2)*"#",string(1/deltacorrect)*"#",string(dsigmadt(sSmol,sBig,tSmol,tBig,uSmol,uBig))*"#")

    if (Sval==Inf || Sval == -Inf || Sval < 0e0)
        println("")
        println("Sval = $Sval")
        println("deltacorrect = $deltacorrect")
        println("dsdt = $(dsigmadt(sSmol,sBig,tSmol,tBig,uSmol,uBig))")
        println("sSmol = $sSmol")
        println("tSmol = $tSmol")
        println("uSmol = $uSmol")
        println("p1v = $p1v")
        println("p2v = $p2v")
        println("p4v = $p4v")
        error("Sval Inf")  
    end

#=     if (Sval/prob) >= 1e5
        println("")
        println("sSmol = $sSmol")
        println("tSmol = $tSmol")
        println("uSmol = $uSmol")
        println("check = $(sSmol+tSmol+uSmol+sBig+tBig+uBig)")
    end =#

    return Sval

end

"""
    InvariantFlux(s,mu12,mu22)

returns the value of the invariant flux with 's' Mandelstram variable and masses 'mass1' and 'mass2'
"""
function InvariantFlux(s::Float64,mu12::Float64,mu22::Float64)

    # sqrt(lambda(s,m1^2,m2^2))/2 = sqrt(s)|p*|
    return sqrt(s^2-2*s*(mu12+mu22)+(mu12-mu22)^2)/2

end

"""
    InvariantFluxSmall(sSmol,mu12,mu22)

returns the value of the invariant flux with smalled 's' Mandelstram variable (sSmol = s - (m1+m2)^2)
"""
function InvariantFluxSmall(sSmol::Float64,m1::Float64,m2::Float64)
    # Better accuracy for small s

    # sqrt(lambda(s,m1^2,m2^2))/2 = sqrt(s)|p*|
    return sqrt(sSmol*(sSmol+4*m1*m2))/2

end

"""
    InvariantFlux2(s,mass12,mass22)

returns the value of the squared invariant flux with 's' Mandelstram variable and masses 'mass1' and 'mass2'
"""
function InvariantFlux2(s::Float64,mu12::Float64,mu22::Float64)

    # lambda(s,m1^2,m2^2)/4 = s|p*|^2
    return (s^2-2*s*(mu12+mu22)+(mu12-mu22)^2)/4

end

"""
    InvariantFluxSmall(sSmol,mass12,mass22)

returns the value of the squared invariant flux with smalled 's' Mandelstram variable (sSmol = s - (m1+m2)^2)
"""
function InvariantFlux2Small(sSmol::Float64,m1::Float64,m2::Float64)
    # Better accuracy for small s

    # lambda(s,m1^2,m2^2)/4 = s|p*|^2
    return (sSmol*(sSmol+4*m1*m2))/4

end


function SValueTest(p3v::Vector{Float64},p4v::Vector{Float64},Tval::Float64,p1v::Vector{Float64},p2v::Vector{Float64},mu1::Float64,mu2::Float64,mu3::Float64,mu4::Float64)

    # define normalise masses
    m1 = mu1
    m2 = mu2 
    m3 = mu3
    m4 = mu4 

    # pre-defining terms for efficiency 
    p1::Float64 = p1v[1]
    p2::Float64 = p2v[1]

    ct1::Float64 = p1v[2] #cospi(p1v[2])
    ct2::Float64 = p2v[2] #cospi(p2v[2]) 

    st1::Float64 = sqrt(1e0-p1v[2]^2) #sinpi(p1v[2])
    st2::Float64 = sqrt(1e0-p2v[2]^2) #sinpi(p2v[2])

    ch1h2::Float64 = cospi(p1v[3]-p2v[3])

    ctheta12::Float64 = ct1*ct2 + ch1h2*st1*st2

    m32 = m3^2
    m42 = m4^2
    m12 = m1^2
    m22 = m2^2
    
    Es1::Float64 = m1 != 0e0 ? (p1^2)/(sqrt(m12+p1^2)+m1) : p1
    Es1s::Float64 = Es1/p1
    E1::Float64 = Es1 + m1
 
    Es2::Float64 = m2 != 0e0 ? (p2^2)/(sqrt(m22+p2^2)+m2) : p2
    Es2s::Float64 = Es2/p2
    E2::Float64 = Es2 + m2

    sBig::Float64 = (m1+m2)^2
    #sSmol::Float64 = 2*(m1*Es2 + m2*Es1 + Es1*Es2 - p1*p2*(ct1*ct2+ch1h2*st1*st2))
    sSmol::Float64 = 2*p1*p2*(-ctheta12 + Es1s*Es2s + m1*Es2s/p1 + m2*Es1s/p2)

    # Sspe anisotropic emission spectrum (to be integrated over d^2p1d^3p3d^3p4). See obsidian note on discrete anisotropic kinetic equation
    val::Float64 = Tval*(sBig+sSmol)/InvariantFluxSmall(sSmol,m3,m4)/pi

    p3::Float64 = p3v[1]
    ct3::Float64 = p3v[2] 
    st3::Float64 = sqrt(1e0-p3v[2]^2)
    ch3h1::Float64 = cospi(p3v[3]-p1v[3])
    ch3h2::Float64 = cospi(p3v[3]-p2v[3])
    Es3::Float64 = m3 != 0e0 ? (p3^2)/(sqrt(m32+p3^2)+m3) : p3
    Es3s::Float64 = Es3/p3

    p4::Float64 = p4v[1]
    ct4::Float64 = p4v[2] 
    st4::Float64 = sqrt(1e0-p4v[2]^2)
    ch4h1::Float64 = cospi(p4v[3]-p1v[3])
    ch4h2::Float64 = cospi(p4v[3]-p2v[3])
    Es4::Float64 = m4 != 0e0 ? (p4^2)/(sqrt(m42+p4^2)+m4) : p4
    Es4s::Float64 = Es4/p4

    theta13 = (ct3*ct1+ch3h1*st3*st1)
    theta23 = (ct3*ct2+ch3h2*st3*st2)
    theta14 = (ct4*ct1+ch4h1*st4*st1)
    theta24 = (ct4*ct2+ch4h2*st4*st2)

    deltacorrect3::Float64 = Es1*p3 - Es3*p1*theta13
    deltacorrect3 += m1*p3 - m3*p1*theta13
    deltacorrect3 += Es2*p3 - Es3*p2*theta23    
    deltacorrect3 += m2*p3 - m3*p2*theta23

    deltacorrect4::Float64 = Es1*p4 - Es4*p1*theta14
    deltacorrect4 += m1*p4 - m4*p1*theta14
    deltacorrect4 += Es2*p4 - Es4*p2*theta24
    deltacorrect4 += m2*p4 - m4*p2*theta24
    # more Float accurate for when p1 and p2 have large order of magnitude difference as sum uses pairwise summation to reduce round of errors
    #sumTerms .= (Es1, -Es3prime*p1*(ct3*ct1+ch3h1*st3*st1), Es2, -Es3prime*p2*(ct3*ct2+ch3h2*st3*st2), m1, -(m3/p3)*p1*(ct3*ct1+ch3h1*st3*st1), m2, -(m3/p3)*p2*(ct3*ct2+ch3h2*st3*st2))
    #deltacorrect = p3*sum_oro(sumTerms)

    #if (stuCheck(sSmol,sBig,tSmol,tBig,uSmol,uBig,m1,m2,m3,m4) == false || tCheck(tSmol,tBig,m1,m2,m3,m4) == false || uCheck(uSmol,uBig,m1,m2,m3,m4) == false)
    #    error("stu check")
    #end

    Sval3 = val*(p3^2/(deltacorrect3*sign(deltacorrect3)))
    Sval4 = val*(p4^2/(deltacorrect4*sign(deltacorrect4)))

    return Sval3,Sval4

end

function SValueTest2(p3v::Vector{Float64},p4v::Vector{Float64},dsigmadt::Function,tSmol::Float64,sSmol::Float64,p1v::Vector{Float64},p2v::Vector{Float64},mu1::Float64,mu2::Float64,mu3::Float64,mu4::Float64)

    # define normalise masses
    m1 = mu1
    m2 = mu2 
    m3 = mu3
    m4 = mu4 

    # pre-defining terms for efficiency 
    p1::Float64 = p1v[1]
    p2::Float64 = p2v[1]

    ct1::Float64 = p1v[2] #cospi(p1v[2])
    ct2::Float64 = p2v[2] #cospi(p2v[2]) 

    st1::Float64 = sqrt(1e0-p1v[2]^2) #sinpi(p1v[2])
    st2::Float64 = sqrt(1e0-p2v[2]^2) #sinpi(p2v[2])

    ch1h2::Float64 = cospi(p1v[3]-p2v[3])

    m32 = m3^2
    m42 = m4^2
    m12 = m1^2
    m22 = m2^2

    Es1::Float64 = m1 != 0e0 ? (p1^2)/(sqrt(m12+p1^2)+m1) : p1
    Es1s::Float64 = Es1/p1
    E1::Float64 = Es1 + m1
 
    Es2::Float64 = m2 != 0e0 ? (p2^2)/(sqrt(m22+p2^2)+m2) : p2
    Es2s::Float64 = Es2/p2
    E2::Float64 = Es2 + m2

    sBig::Float64 = (m1+m2)^2
    #sSmol::Float64 = 2*(m1*Es2 + m2*Es1 + Es1*Es2 - p1*p2*(ct1*ct2+ch1h2*st1*st2))
    #sSmol::Float64 = 2*p1*p2*(-ct1*ct2 -ch1h2*st1*st2 + Es1s*Es2s + m1*Es2s/p1 + m2*Es1s/p2)

    # Sspe anisotropic emission spectrum (to be integrated over d^2p1d^3p3d^3p4). See obsidian note on discrete anisotropic kinetic equation
    #val::Float64 = (1/E1)*(1/E2)*(sBig+sSmol)*InvariantFluxSmall(sSmol,m1,m2)/InvariantFluxSmall(sSmol,m3,m4)/(pi)
    val::Float64 = (1/E1)*(1/E2)*InvariantFlux2Small(sSmol,m1,m2)/(pi)

    p3::Float64 = p3v[1]
    t3::Float64 = p3v[4]
    (st3::Float64,ct3::Float64) = sincospi(t3)
    #ct3::Float64 = p3v[2] 
    #st3::Float64 = sqrt(1e0-p3v[2]^2)
    ch3h1::Float64 = cospi(p3v[3]-p1v[3])
    ch3h2::Float64 = cospi(p3v[3]-p2v[3])
    Es3::Float64 = m3 != 0e0 ? (p3^2)/(sqrt(m32+p3^2)+m3) : p3
    Es3s::Float64 = Es3/p3

    theta13 = (ct3*ct1+ch3h1*st3*st1)
    theta23 = (ct3*ct2+ch3h2*st3*st2)

    # t = tBig + tSmol
    tBig3::Float64 = (m3-m1)^2
    #tSmol3::Float64 = 2*p3*p1*(theta13 - Es3s*Es1s - m1*Es3s/p1 - m3*Es1s/p3)
    # u = uBig + uSmol
    uBig3::Float64 = (m2-m3)^2
    #uSmol3::Float64 = -2*p2*p3*(-ct2*ct3 -ch3h2*st2*st3 + Es2s*Es3s + m3*Es2s/p3 + m2*Es3s/p2)
    #println("uSmol3 = $uSmol3")

    #uSmol3 = m1^2+m2^2+m3^2+m4^2-sBig-tBig3-uBig3-sSmol-tSmol3 
    #println("uSmol3 = $uSmol3") 
    uSmol3 = m1^2+m2^2+m3^2+m4^2-sBig-tBig3-uBig3-sSmol-tSmol 
    #println("uSmol3 = $uSmol3") 

    deltacorrect3::Float64 = Es1*p3 - Es3*p1*theta13
    deltacorrect3 += m1*p3 - m3*p1*theta13
    deltacorrect3 += Es2*p3 - Es3*p2*theta23    
    deltacorrect3 += m2*p3 - m3*p2*theta23

    Sval3::Float64 = dsigmadt(sSmol,sBig,tSmol,tBig3,uSmol3,uBig3)*val*(p3^2/(deltacorrect3*sign(deltacorrect3)))

    p4::Float64 = p4v[1]
    ct4::Float64 = p4v[2] 
    st4::Float64 = sqrt(1e0-p4v[2]^2)
    ch4h1::Float64 = cospi(p4v[3]-p1v[3])
    ch4h2::Float64 = cospi(p4v[3]-p2v[3])
    Es4::Float64 = m4 != 0e0 ? (p4^2)/(sqrt(m42+p4^2)+m4) : p4
    Es4s::Float64 = Es4/p4

    # t = tBig + tSmol
    tBig4::Float64 = (m2-m4)^2
    #tSmol4::Float64 = 2*p2*p4*(ct2*ct4 +ch4h2*st2*st4 - Es2s*Es4s - m4*Es2s/p4 - m2*Es4s/p2)
    # u = uBig + uSmol
    uBig4::Float64 = (m4-m1)^2
    #uSmol::Float64 = -2*(m1*Es4 + m4*Es1 + Es4*Es1 - p4*p1*(ct4*ct1+ch4h1*st4*st1))
    #uSmol4::Float64 = -2*p4*p1*(-ct4*ct1 -ch4h1*st4*st1 + Es4s*Es1s + m1*Es4s/p1 + m4*Es1s/p4)

    uSmol4 = m1^2+m2^2+m3^2+m4^2-sBig-tBig3-uBig3-sSmol-tSmol 

    theta14 = (ct4*ct1+ch4h1*st4*st1)
    theta24 = (ct4*ct2+ch4h2*st4*st2)

    deltacorrect4::Float64 = Es1*p4 - Es4*p1*theta14
    deltacorrect4 += m1*p4 - m4*p1*theta14
    deltacorrect4 += Es2*p4 - Es4*p2*theta24
    deltacorrect4 += m2*p4 - m4*p2*theta24

    Sval4::Float64 = dsigmadt(sSmol,sBig,tSmol,tBig4,uSmol4,uBig4)*val*(p4^2/(deltacorrect4*sign(deltacorrect4)))

    #=t_error = 1-tSmol3/tSmol4
    if t_error > 1e2 || t_error < -1e2
        println("t_error: $t_error")
        println("ts3=$tSmol3")
        println("ts4=$tSmol4")
        println("us3=$uSmol3")
        println("us4=$uSmol4")
        println("p1=$p1")
        println("p2=$p2")
        println("p3=$p3")
        println("p4=$p4")
        println("theta13 = $theta13")
        println("sct1 = $st1,$ct1")
        println("sct3 = $st3,$ct3")
        println("ch3h1 = $ch3h1")
        println("theta24 = $theta24")
        println("st2,ct2 = $st2,$ct2")
        println("st4,ct4 = $st4,$ct4")
        println("ch4h2 = $ch4h2")
        println("")
    end=#

    # more Float accurate for when p1 and p2 have large order of magnitude difference as sum uses pairwise summation to reduce round of errors
    #sumTerms .= (Es1, -Es3prime*p1*(ct3*ct1+ch3h1*st3*st1), Es2, -Es3prime*p2*(ct3*ct2+ch3h2*st3*st2), m1, -(m3/p3)*p1*(ct3*ct1+ch3h1*st3*st1), m2, -(m3/p3)*p2*(ct3*ct2+ch3h2*st3*st2))
    #deltacorrect = p3*sum_oro(sumTerms)

    #if (stuCheck(sSmol,sBig,tSmol,tBig,uSmol,uBig,m1,m2,m3,m4) == false || tCheck(tSmol,tBig,m1,m2,m3,m4) == false || uCheck(uSmol,uBig,m1,m2,m3,m4) == false)
    #    error("stu check")
    #end

    return Sval3,Sval4

end

#=
tS = -0.2844979694423338
t_error : -408.8789407386356
ts3 = -108.9888855988345
ts4 = -0.2659050630960145
us3 = -1.1204686725550848
us4 = -0.5155616751268332
p1 = 10.890835664944595
p2 = 0.04332837367344047
p3 = 6.915797711629199
p4 = 3.992255891330024
theta13 = 0.27785649953134167
st1,ct1 = 0.9689484526950011, 0.24726280759540686
st3,ct3 = 0.03170243009365831, 0.9994973516353891
ch3h1 = 0.999998619568608
m3 = 1.0
m32= m3^2
m1 = 1.0
m12 = m1^2

Es1::Float64 = m1 != 0e0 ? (p1^2)/(sqrt(m12+p1^2)+m1) : p1
Es1s::Float64 = Es1/p1

Es3::Float64 = m3 != 0e0 ? (p3^2)/(sqrt(m32+p3^2)+m3) : p3
Es3s::Float64 = Es3/p3

theta13 = (ct3*ct1+ch3h1*st3*st1)
theta23 = (ct3*ct2+ch3h2*st3*st2)

# t = tBig + tSmol
tBig3::Float64 = (m3-m1)^2
tSmol3::Float64 = 2*p3*p1*(theta13 - Es3s*Es1s - m1*Es3s/p1 - m3*Es1s/p3)

tSmoltest = 2*m1*m3-2*(Es1+m1)*(Es3+m3) +2*p1*p3*theta13
=#

#=
sSmol = 0.0021691691589844114
sBig = 1.0
tSmol = -4.8624970987408064e-6
tBig = 0.0
uSmol = -0.002164229519333684
uBig = 1.0
s+t+u = 2.000000077142552
m1 = 1.0
m2 = 0.0
m3 = 1.0
m4 = 0.0
p1v = [38316.52012320547, 0.45880103655549265, 1.8997882467761282]
p2v = [5.315929608854946e-8, 0.23648605336493644, 0.2632690527850132]
p3v = [38229.26560868081, 0.4588010365570692, 1.8997882467772347, 0.3482790688409888]

    # pre-defining terms for efficiency 
    p1::Float64 = p1v[1]
    p2::Float64 = p2v[1]

    ct1::Float64 = p1v[2] #cospi(p1v[2])
    ct2::Float64 = p2v[2] #cospi(p2v[2]) 

    st1::Float64 = sqrt(1e0-p1v[2]^2) #sinpi(p1v[2])
    st2::Float64 = sqrt(1e0-p2v[2]^2) #sinpi(p2v[2])

    ch1h2::Float64 = cospi(p1v[3]-p2v[3])

    m32 = m3^2
    m42 = m4^2
    m12 = m1^2
    m22 = m2^2
    
    Es1::Float64 = m1 != 0e0 ? (p1^2)/(sqrt(m12+p1^2)+m1) : p1
    Es1s::Float64 = Es1/p1
    E1::Float64 = Es1 + m1
 
    Es2::Float64 = m2 != 0e0 ? (p2^2)/(sqrt(m22+p2^2)+m2) : p2
    Es2s::Float64 = Es2/p2
    E2::Float64 = Es2 + m2

    sBig = (m1+m2)^2
    #sSmol::Float64 = 2*(m1*Es2 + m2*Es1 + Es1*Es2 - p1*p2*(ct1*ct2+ch1h2*st1*st2))
    sSmol = 2*p1*p2*(-ct1*ct2 + Es1s*Es2s + m1*Es2s/p1 + m2*Es1s/p2 - ch1h2*st1*st2)

    # Sspe anisotropic emission spectrum (to be integrated over d^2p1d^3p3d^3p4). See obsidian note on discrete anisotropic kinetic equation
    val::Float64 = (1/E1)*(1/E2)*(InvariantFlux2Small(sSmol,m1,m2))/pi

    p3::Float64 = p3v[1]
    ct3::Float64 = p3v[2] # sinpi and cospi slightly slower than sin(pi*) but more accurate apparently
    st3::Float64 = sqrt(1e0-p3v[2]^2)
    ch3h1::Float64 = cospi(p3v[3]-p1v[3])
    ch3h2::Float64 = cospi(p3v[3]-p2v[3])
    Es3::Float64 = m3 != 0e0 ? (p3^2)/(sqrt(m32+p3^2)+m3) : p3
    Es3s::Float64 = Es3/p3

    # t = tBig + tSmol
    tBig = (m3-m1)^2
    #tSmol::Float64 = -2*(m1*Es3 + m3*Es1 + Es3*Es1 - p3*p1*(ct3*ct1+ch3h1*st3*st1))
    tSmol = 2*p3*p1*(ct3*ct1 - Es3s*Es1s - m1*Es3s/p1 - m3*Es1s/p3 + ch3h1*st3*st1)
    # u = uBig + uSmol
    uBig = (m2-m3)^2
    #uSmol::Float64 = m12+m22+m32+m42 - sBig - tBig - uBig - sSmol - tSmol  # this leads to Float64 issues better to calculate directly
    #uSmol::Float64 = -2*(m3*Es2 + m2*Es3 + Es2*Es3 - p2*p3*(ct2*ct3+ch3h2*st2*st3))
    uSmol = 2*p2*p3*(ct2*ct3 - Es2s*Es3s - m3*Es2s/p3 - m2*Es3s/p2 + ch3h2*st2*st3)

    deltacorrect::Float64 = (Es1*p3 - Es3*p1*(ct3*ct1+ch3h1*st3*st1) - m3*p1*(ct3*ct1+ch3h1*st3*st1)) + m1*p3 + m2*p3 + (Es2*p3 - Es3*p2*(ct3*ct2+ch3h2*st3*st2) - m3*p2*(ct3*ct2+ch3h2*st3*st2))
    
    (-(Es2+m2)-(ct3*ct2+ch3h2*st3*st2)*p2+m1^2/2/p1)*p3-m3^3*p1/p3/2-(ct3*ct2+ch3h2*st3*st2)*m3^2*p2/p3/2


pa = (p1v[1]+p3v[1])/2
pd = (p1v[1]-p3v[1])/2
m1 = 1e0
m3 = 1e0

-2(pa^2 * (-1 + p1v[2]*p3v[2] + cospi(p1v[3]-p3v[3])*sqrt(1-p1v[2]^2)*sqrt(1-p3v[2]^2)))

+2(-1 + (2*m3^2)/(m3^2 + pa^2) + p1v[2]*p3v[2] + cospi(p1v[3]-p3v[3])*sqrt(1-p1v[2]^2)*sqrt(1-p3v[2]^2))*pd^2

p1v[2]*p3v[2] + cospi(p1v[3]-p3v[3])*sqrt(1-p1v[2]^2)*sqrt(1-p3v[2]^2)

a = (acos(p1v[2])+acos(p3v[2]))/2
b = (acos(p1v[2])-acos(p3v[2]))/2
c1 = pi*(p1v[3]-p3v[3])/2
b^2*(-2 + 2c1^2) - 2c1^2*sin(a)^2

acos(p1v[2])-pi*p3v[4]

test  = -2((p1v[1]+p3v[1])^2/4 * (b^2*(-2 + 2c1^2) - 2c1^2*sin(a)^2))
test += +2((2*m3^2)/(m3^2 + pa^2) + b^2*(-2 + 2c1^2) - 2c1^2*sin(a)^2)*pd^2
sSmol+uSmol-test

=#

p3v = [8.492263520684787, -0.9728771747352043, 0.9919779470075307, 0.9256948160192373]
p1v = [8.492263998166756, -0.9728771750199965, 0.9919779468286327, 0.9256948164111243]
p2v = [7.463234147568588e-10, -0.9744986069019634, 0.7018819330469559, 0.9279598929075087]

   # pv should be [p,t,h]
   p1::Float64 = p1v[1]
   p2::Float64 = p2v[1]

   st1::Float64,ct1::Float64 = sincospi(p1v[4])
   st2::Float64,ct2::Float64 = sincospi(p2v[4])
   st3::Float64,ct3::Float64 = sincospi(p3v[4])

   ch3h1::Float64 = cospi(p3v[3]-p1v[3])
   ch3h2::Float64 = cospi(p3v[3]-p2v[3])
   ch1h2::Float64 = cospi(p1v[3]-p2v[3])

   m3 = 1.0
    m4 = 0.0
    m1 = 1.0
    m2 = 0.0
   m32::Float64 = m3^2
   m42::Float64 = m4^2
   m12::Float64 = m1^2
   m22::Float64 = m2^2

   p12::Float64 = p1^2
   p22::Float64 = p2^2

   E1::Float64 = sqrt(m12+p12)
   E2::Float64 = sqrt(m22+p22)

   ctheta12::Float64 = ct1*ct2 + ch1h2*st1*st2
   ctheta13::Float64 = ct3*ct1 + ch3h1*st3*st1
   ctheta23::Float64 = ct3*ct2 + ch3h2*st3*st2

   C3sqr::Float64 = (E1+E2)^2*((m32-m42+m12+m22+2*E1*E2-2*p1*p2*ctheta12)^2-4*m32*((E1+E2)^2-(p1*ctheta13+p2*ctheta23)^2))

   C2::Float64 =-4*(p1*ctheta13+p2*ctheta23)*(m32-m42+m12+m22+2*E1*E2-2*p1*p2*ctheta12)

   C4::Float64 = -8*((E1+E2)^2-(p1*ctheta13+p2*ctheta23)^2)

   C3 = 4*sqrt(C3sqr)

   (C2-C3)/C4
   (C2+C3)/C4

   E1 = sqrt(m12+p12)
   E2 = sqrt(m22+p22)
   p3 = p3v[1]
   E3 = sqrt(m32+p3v[1]^2)

   Es1::Float64 = m1 != 0e0 ? (p1^2)/(sqrt(m12+p1^2)+m1) : p1
   Es1s::Float64 = Es1/p1
   E1::Float64 = Es1 + m1

   Es2::Float64 = m2 != 0e0 ? (p2^2)/(sqrt(m22+p2^2)+m2) : p2
   Es2s::Float64 = Es2/p2
   E2::Float64 = Es2 + m2

   Es3::Float64 = m3 != 0e0 ? (p3^2)/(sqrt(m32+p3^2)+m3) : p3
   Es3s::Float64 = Es3/p3
   E3 = Es3 + m3

   tSmol::Float64 = 2*p3*p1*(ctheta13 - Es3s*Es1s - m1*Es3s/p1 - m3*Es1s/p3)

   deltacorrect::Float64 = Es1*p3 - Es3*p1*ctheta13
    deltacorrect += m1*p3 - m3*p1*ctheta13
    deltacorrect += Es2*p3 - Es3*p2*ctheta23    
    deltacorrect += m2*p3 - m3*p2*ctheta23

    E1*p3 - E3*p1*ctheta13 + E2*p3 - E3*p2*ctheta23    

