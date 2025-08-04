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
    LossValue(p1v,p2v,sigma,mu1,mu2)

returns `LossVal` based on initial momentum states `p1v` and `p2v` and cross section `sigma` based on particle selection.
```math
Loss_\\text{val} = \\frac{1}{p^0_1p^0_2}\\sigma(s)F_12(s)
```
If initial state fails `sCheck`, i.e. cannot generate a physical output state, LossVal is set to 0e0. 
Assumes f(x,p,u,ϕ)=f(x,vec{p})p^2=constant over bin
"""
function LossValue(p1v::Vector{Float64},p2v::Vector{Float64},sigma::Function,m1::Float64,m2::Float64,m3::Float64,m4::Float64)

    # pre-defining terms for efficiency 
    p1::Float64 = p1v[1]
    p2::Float64 = p2v[1]

    st1::Float64,ct1::Float64 = sincospi(p1v[4])
    st2::Float64,ct2::Float64 = sincospi(p2v[4])

    ch1h2::Float64 = cospi(p1v[3]-p2v[3])

    ctheta12::Float64 = ct1*ct2+ch1h2*st1*st2

    m12::Float64 = m1^2
    m22::Float64 = m2^2

    LossVal::Float64 = 0e0
    
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
        
        LossVal = (1/E1)*(1/E2)*(InvariantFluxSmall(sSmol,m1,m2))*sigma(sSmol,sBig)
        if LossVal==Inf || LossVal < 0e0
            println("")
            println("p1v = $p1v")
            println("p2v = $p2v")
            println("LossVal = $LossVal")
            println("sSmol = $sSmol")
            println("sBig = $sBig")
            error("ST1 Inf or -ve#")
        end
    else # if not valid set T value to zero
        LossVal = 0e0
    end

    return LossVal, sBig, sSmol

end


"""
    GainValue3(p3v,p1v,p2v,dsigmadt,mu1,mu2,mu3)

Returns `GainVal` based on initial momentum states `p1v` and `p2v` and final state `p3v` and differential cross section `dsigmadt` based on particle selection 12->34. 
```math
Gain_\\text{val}=\\frac{\\mathrm{d}\\sigma_{12|34}}{\\mathrm{d}t}\\frac{\\mathcal{F}_{12}^2}{\\pi\\left|p^0_3(p_1\\cos\\Theta_{31}+p_2\\cos\\Theta_{32})-p_3(p^0_1+p^0_2)\\right|}\frac{p_3^2}{p^0_1p^0_2}.
``` 
Assumes f(x,p,u,ϕ)=f(x,vec{p})p^2=constant over bin
"""
function GainValue3(p3v::Vector{Float64},p1v::Vector{Float64},p2v::Vector{Float64},sBig::Float64,sSmol::Float64,dsigmadt::Function,m1::Float64,m2::Float64,m3::Float64,m4::Float64)

    # pre-defining terms for efficiency 
    p1::Float64 = p1v[1]
    p2::Float64 = p2v[1]
    st1::Float64,ct1::Float64 = sincospi(p1v[4])
    st2::Float64,ct2::Float64 = sincospi(p2v[4])

    ch1h2::Float64 = cospi(p1v[3]-p2v[3])

    p3::Float64 = p3v[1]
    st3::Float64,ct3::Float64 = sincospi(p3v[4])
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

    val::Float64 = (1/E1)*(1/E2)*(InvariantFlux2Small(sSmol,m1,m2))/pi

    #=if (stuCheck(sSmol,sBig,tSmol,tBig,uSmol,uBig,m1,m2,m3,m4) == false || tCheck(tSmol,tBig,m1,m2,m3,m4) == false || uCheck(uSmol,uBig,m1,m2,m3,m4) == false)
        println("p1v = $p1v")
        println("p2v = $p2v")
        println("p3v = $p3v")
        error("stu check")
    end=#

    GainVal = dsigmadt(sSmol,sBig,tSmol,tBig,uSmol,uBig)*val*(p3^2/(deltacorrect*sign(deltacorrect)))

    if (GainVal==Inf || GainVal == -Inf || GainVal < 0e0)
        println("")
        println("GainVal = $GainVal")
        println("deltacorrect = $deltacorrect")
        println("dsdt = $(dsigmadt(sSmol,sBig,tSmol,tBig,uSmol,uBig))")
        println("sSmol = $sSmol")
        println("tSmol = $tSmol")
        println("uSmol = $uSmol")
        println("p1v = $p1v")
        println("p2v = $p2v")
        println("p3v = $p3v")
        error("GainVal Inf")  
    end

    return GainVal

end


"""
    GainValue4(p3v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)

Returns `GainVal` based on initial momentum states `p1v` and `p2v` and final state `p4v` and differential cross section `dsigmadt` based on particle selection 12->34.  
Assumes f(x,p,μ)=constant over bin
"""
function GainValue4(p4v::Vector{Float64},p1v::Vector{Float64},p2v::Vector{Float64},sBig::Float64,sSmol::Float64,dsigmadt::Function,m1::Float64,m2::Float64,m3::Float64,m4::Float64,prob::Float64)

    # pre-defining terms for efficiency 
    p1::Float64 = p1v[1]
    p2::Float64 = p2v[1]
    st1::Float64,ct1::Float64 = sincospi(p1v[4])
    st2::Float64,ct2::Float64 = sincospi(p2v[4])

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

    # u = uBig + uSmol
    uBig::Float64 = (m4-m1)^2
    #uSmol::Float64 = -2*(m1*Es4 + m4*Es1 + Es4*Es1 - p4*p1*(ct4*ct1+ch4h1*st4*st1))
    uSmol::Float64 = 2*p4*p1*(ctheta14 - Es4s*Es1s - m1*Es4s/p1 - m4*Es1s/p4)
    # t = tBig + tSmol
    tBig::Float64 = (m2-m4)^2
    #tSmol::Float64 = m12+m22+m32+m42 - sBig - uBig - tBig - sSmol - uSmol # Leads to Float64 issues, better to calculate directly 
    #tSmol::Float64 = -2*(m4*Es2 + m2*Es4 + Es2*Es4 - p2*p4*(ct2*ct4+ch4h2*st2*st4))
    tSmol::Float64 = 2*p2*p4*(ctheta24 - Es2s*Es4s - m4*Es2s/p4 - m2*Es4s/p2)
    
    val::Float64 = (1/E1)*(1/E2)*(InvariantFlux2Small(sSmol,m1,m2))/pi

    #=if (stuCheck(sSmol,sBig,tSmol,tBig,uSmol,uBig,m1,m2,m3,m4) == false || tCheck(tSmol,tBig,m1,m2,m3,m4) == false || uCheck(uSmol,uBig,m1,m2,m3,m4) == false)
        error("stu check")
    end=#

    GainVal = dsigmadt(sSmol,sBig,tSmol,tBig,uSmol,uBig)*val*(p4^2/(deltacorrect*sign(deltacorrect)))

    if (GainVal==Inf || GainVal == -Inf || GainVal < 0e0)
        println("")
        println("GainVal = $GainVal")
        println("deltacorrect = $deltacorrect")
        println("dsdt = $(dsigmadt(sSmol,sBig,tSmol,tBig,uSmol,uBig))")
        println("sSmol = $sSmol")
        println("tSmol = $tSmol")
        println("uSmol = $uSmol")
        println("p1v = $p1v")
        println("p2v = $p2v")
        println("p4v = $p4v")
        error("GainVal Inf")  
    end

    return GainVal

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


  

