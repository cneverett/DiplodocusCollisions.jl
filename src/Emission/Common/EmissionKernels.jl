"""
    SyncKernel(p3v,p1v,m1,z1,B)

Returns the emission rate for a single photon ``p3v`` state emitted by a charged particle in state ``p1v`` with charge ``z1`` relative to the fundamental charge and mass ``m1`` relative to the mass of the electron, in a uniform magnetic field ``B``.
"""
function SyncKernel(p3v::Vector{Float64},p1v::Vector{Float64},m1::Float64,z1::Float64,Ext::Vector{Float64})

    # p3is Photon
    # p1 is Charged Particl
    
    B::Float64 = Ext[1] # B field in Tesla

    p3::Float64 = p3v[1]
    p1::Float64 = p1v[1]
    st1::Float64,ct1::Float64 = sincospi(p1v[4])
    st3::Float64,ct3::Float64 = sincospi(p3v[4])
    
    E1::Float64 = sqrt(p1^2 + m1^2)

    Jfactor1 = (E1*ct3-p1*ct1)/(st3) # code breaks if st3 = 0 FIX
    Jfactor2 = p1*st1

    n = abs((mEle^2*c^2)/(z1*ħ*q*B)) * p3 * (E1-p1*ct3*ct1)
    #println(n)

    y = p1 * st1 *st3 / (E1-p1*ct1*ct3) # y=x/n
    #println(y)

    # characteristic frequency
    ω0 = abs((z1*q*B))/(E1*mEle)
    #println("critical photon momentum: "*string(ħ*ω0/(mEle*c^2)*E1^3))


    if n > 1e2 #&& 1-y < 0.01
        # approximation for J's to second order 
        e = 1-y^2
        K13 = besselk(1/3,n*e^(3/2)/3)
        K23 = besselk(2/3,n*e^(3/2)/3)
        J1 = ((sqrt(e))/(pi*sqrt(3)))*(K13 #=+(e/10)*(K13-2*n*e^(3/2)*K23)=#)
        #println(K23)
        J2 = (e/(pi*sqrt(3)))*(K23 #=+ (e/5)*(2*K23-(1/(e^(3/2)*n)+n*e^(3/2))*K13)=#)
    elseif n < 1e0
        # omega < omega0 therefore no synchrotron radiation
        J1 = 0e0
        J2 = 0e0
    else
        # exact J's
        J1 = besselj(n,n*y)
        J2 = 1/2 * (besselj(n-1,n*y) - besselj(n+1,n*y))
    end

    val = (p3/E1)*((Jfactor1*J1)^2+(Jfactor2*J2)^2)
    #println(val)

    factor = (abs(z1/B))*(3*c^4*mEle^5)/(4*pi*ħ^3*μ0*q^3) # synchrotron emission rate divided by c*σT

    #println(factor)

    return val*factor
    
end

#SyncKernel([1e-14,0.6],[1e1,0.5],1,1,1e-4)
#n = 7.36e10
#e = 1-(0.0045)^2
#besselk(1/3,n*e^(3/2)/3)