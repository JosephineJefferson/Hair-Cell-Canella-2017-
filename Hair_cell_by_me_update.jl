# Hair cell model

# Based on publication: Canella et al. (2017)

# Josephine Jefferson
# May 2020

#Currents:
# IKD = compound delayed K current
#(above comprises v-gated K (IKV), Ca gated K (IKCa), and Ca (ICa))
#IKV and IKCa do not have their own equations and are functionally indiscernable
#ICa is important for neurotransmitter release, and can be modelled separately

#IA = fast transient K current
#IL = leak current
#I_trans = transduction current (mechanoceptors)
using Unitful
"""
numbers here for reference

const g_KD_max = 0.01u"μS" #maximal KD conductance
const g_A_max = 0.02u"μS" #maximal A (K-based) conductance
const g_Ca_max = 1.2u"μS" #maximal Ca conductance
const P_Ca = 1.01e-10u"cm/s" #Ca permeability
const g_L = 0.77u"nS" #leakage conductance (constant)

#g_trans (transduction conductance) is described on the paper (how to calc)

const E_K = -96.0u"mV" #Postassium equilibrium Potential
const E_Ca = 40.0u"mV" #Apparent ICa equilibrium Potential
const E_L = -68.9u"mV" #Leak equilibrium Potential
const E_trans = 0.0u"mV" #transduction equilibrium Potential

#tau_C and tau_Ca (time constants for ICA activation - sinusoidal and step
#model) are defined in paper

# time constants of ICA decay
const tau_2 = 600u"ms"
#const tau_1 = 6000*exp(-t/tau_2) #dependent on time
"""
#CONSTANT TIMESTEP
const dt=1e-11u"s"

#will probably want to 1/0-proof these
#V-dependent parameters of IKD current (delayed K and Ca current)
a_KD(V)=1/(1+exp((-14.9-(ustrip(V)))/16.7))
tau_KD(V)=0.89+0.96*exp(-0.5((ustrip(V)-30)/17.3)^2)

#V-dependent parameters of IA current (fast K current)
a_A(V)= 1/(1+exp((-5.7-ustrip(V))/22.2))
h_A(V)= 1/(1+exp((-70.2-ustrip(V))/9.2))
tau_aA(V)= 0.66+1.85*exp(-0.5(abs(ustrip(V)-39.9)/19.4)^1.46)
tau_hA(V) = 22.85*exp(-exp(-(ustrip(V)+39.59)/17.2)-(ustrip(V)+39.59)/17.2 + 1) + 18.4

#a_inf and h_inf
#kluged the constants bc honestly idk and they are not in the paper
const V_a=-68.0u"mV" #chose
const K_a=2.0 #honestly random
const V_h=-10.0u"mV" #chose
const K_h=2.0 #honestly random
a_inf(V)=1/(1+exp((ustrip(V_a)-ustrip(V))/K_a)) #true
h_inf(V)=1/(1+exp((ustrip(V_h)-ustrip(V))/K_h))

#Calculate change of h and a values
Δa_KD(V)=(a_inf(V)-a_KD(V))*(1-exp(-ustrip(dt)/tau_KD(V)))
Δa_A(V)=(a_inf(V)-a_A(V))*(1-exp(-ustrip(dt)/tau_aA(V)))
Δh_A(V)=(h_inf(V)-h_A(V))*(1-exp(-ustrip(dt)/tau_hA(V)))

#Calculate new values based on old values and change
a_KD_new(V, pre)=ustrip(pre)+Δa_KD(V)
a_A_new(V, pre)=ustrip(pre)+Δa_A(V)
h_A_new(V, pre)=ustrip(pre)+Δh_A(V)

"""
So far I am not isolating Ca2+ conductance, and I am not modelling
transduction. IE I am just looking at membrane properties with regards
to voltage (Ca is lumped in with other positive, will want to bring out
later due to relevance to transmitter release).
"""

struct Hair_Cell
    #state vector
    x::Array{Any,1} #changing values #V, IKD, IA, IL, #I_trans(later)
    mem::Array{Any,1} #holds previous a_KD, a_A, h_A
    #Parameters (fixed once declared)
    #Conductances (+poss Ca permeability)
    A::typeof(1.0u"cm^2") #area
    C::typeof(1.0u"μF") #capacitance
    g_KD_max::typeof(1.0u"μS") #maximal KD conductance
    g_A_max::typeof(1.0u"μS") #maximal A (K-based) conductance
    #g_Ca_max::typeof(1.0u"μS") #maximal Ca conductance
    #P_Ca::typeof(1.0u"cm/s") #Ca permeability
    g_L::typeof(1.0u"nS") #leakage conductance (constant)
    #Equilibrium potentials
    E_K::typeof(1.0u"mV") #Postassium equilibrium Potential
    #E_Ca::typeof(1.0u"mV") #Apparent ICa equilibrium Potential
    E_L::typeof(1.0u"mV") #Leak equilibrium Potential
    E_trans::typeof(1.0u"mV") #transduction equilibrium Potential
end

function Hair_Cell(V::typeof(1.0u"mV"), d::typeof(1.0u"nm"))
    Cs=1.0u"μF/cm^2" #specific capacitance(of lipid bilayer)
    #Based on neuron as I could not find a value in Canella paper
    #Constant Parameters
    #Conductances (+ poss Ca permeability)
    A=uconvert(u"cm^2",π*(d^2)) #membrane area calculated based on diameter (passed to constructor) - in cm^2
    C=uconvert(u"μF",Cs*A) #c based on neuron Cs and hc size
    g_KD_max = 0.01u"μS" #maximal KD conductance (delayed K and Ca)
    g_A_max = 0.02u"μS" #maximal A (fast K-based) conductance
    #g_Ca_max = 1.2u"μS" #maximal Ca conductance
    #P_Ca = 1.01e-10u"cm/s" #Ca permeability
    g_L = 0.77u"nS" #leakage conductance (constant)
    #Equilibrium potentials
    E_K = -96.0u"mV" #Postassium equilibrium Potential
    #E_Ca = 40.0u"mV" #Apparent ICa equilibrium Potential
    E_L = -68.9u"mV" #Leak equilibrium Potential
    E_trans = 0.0u"mV" #transduction equilibrium Potential

    #Initial values for channel currents (currently setting to zero rather than calculating)
    IKD = 0.0u"nA" #uconvert(u"nA",g_KD_max*(V-E_K)*a_KD(V))#*tau_KD(V))
    IA =  0.0u"nA" #uconvert(u"nA",g_A_max*(V-E_K)*a_A(V)*h_A(V))#*tau_aA(V)*tau_hA(V))
    IL =  0.0u"nA" #uconvert(u"nA",g_L*(V-E_L))

    #Intial values for a_KD, a_A, h_A
    a_KD=a_inf(V)
    a_A=a_inf(V)
    h_A=h_inf(V)

    return Hair_Cell([V,IKD,IA,IL],[a_KD,a_A,h_A], A,C,g_KD_max,g_A_max,g_L,E_K,E_L,E_trans)
end

"""
Function taking a hair cell and input current and timestep as input.
Alters state vector of hair cell based on input
"""
function update(cell::Hair_Cell,input::typeof(1.0u"nA"))
    v_old=cell.x[1]
    IKD=cell.x[2]
    IA=cell.x[3]
    IL=cell.x[4]

    V= uconvert(u"mV", v_old - dt*((IKD+IA+IL)-input)/cell.C)

    #DO NOT NEED, but this is what these things are
    #a_KD_old=cell.mem[1]
    #a_A_old=cell.mem[2]
    #h_A_old=cell.mem[3]

    #update state
    cell.mem[1]=a_KD_new(V,cell.mem[1])
    cell.mem[2]=a_A_new(V,cell.mem[2])
    cell.mem[3]=h_A_new(V,cell.mem[3])
    cell.x[1]=V
    cell.x[2] = uconvert(u"nA",cell.g_KD_max*(V-cell.E_K)*cell.mem[1])
    cell.x[3] = uconvert(u"nA",cell.g_A_max*(V-cell.E_K)*cell.mem[2]*cell.mem[3])
    cell.x[4] = uconvert(u"nA",cell.g_L*(V-cell.E_L))
end


bella=Hair_Cell(-10.0u"mV",200.0u"nm") #Lysakowski and goldberg 1997: chinchilla hair cell diam 90-400nm, cited in Eatock and Fay
println(bella.x)
for i in 1:1000
    update(bella,0.0u"nA")
    println(bella.x[1])
end
println(bella.x)
