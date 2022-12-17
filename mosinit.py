
from scipy.stats import norm
from neuron import h, load_mechanisms
from numpy import trapz
cvode = h.CVode()
cvode.active(1)
cvode.maxstep(0.2)
h.load_file('stdlib.hoc')
import pylab
pylab.interactive(1) 
import matplotlib.cm as cm
celsius = 36.0 # temperature
Epas = -69

## SELECT MORPHOLOGY
h("forall delete_section()")
h("sec_counted=0")
h.load_file(1,"Morf_default.hoc") # Default, from Halnes2011
h.load_file(1,"fixnseg.hoc") # Segmentize (technicality)

#################################################################
## CaT DISTRIBUTION. (Tdist)
dist = {0:'Null', 1: 'Soma', 2 : 'Proximal', 3 : 'Uniform', 4: 'Middle', 5 : 'Linear',  6:'Distal', 7: 'Halnes'} 
colour = {'Null':'k','Soma':'r', 'Proximal':'y','Uniform': 'b', 'Middle':'c', 'Linear': 'g', 'Distal':'m'}

Tdist = 3 # Choose distribution number between 0 and 7

##################################################################
### Current stimulus
# Choose type of signal input:
# 1 for 1000 ms weak input signal; 2 for 10 ms strong pulse; 3 for synaptic input
input_cond = 1 # Must be 1, 2 or 3

# den is only used when input_cond == 3
den = 99 # Section at which synaptic input is injected (distalmost dendrite): 

if (input_cond == 1):
    dur = 1000
    istim = 0.036
    gna = 0.18
    dsynweight = 0 #0.07341 #0.001 #  Distal synapses, weight (muS)

elif (input_cond == 2):
    dur = 10
    istim = 0.25
    gna = 0.18
    dsynweight = 0 #0.07341 #0.001 #  Distal synapses, weight (muS)
     
elif (input_cond == 3):
    dur = 0
    istim = 0
    gna = 0.18
    dsynweight = 0.002 #0.07341 #0.001 #  Distal synapses, weight (muS)   

else:
    print("Input condition must be 1, 2 or 3")
                  
######################################################################
def test(Tdist,dur,gna,istim,dsynweight):
    global den
    rall = 113 # most people find something like this
    cap = 1.1
    Rm = 45000.0
    Vrest = - 69

    # Channel densities & shifts
    # gna = 0.15
    nash = - 50.5 #-51.3 #-50.3 
    gkdr = 0.4
    kdrsh = - 49.3
    gahp = 0 #0.00013
    gcat = 8.5e-6
    gcal = 0 #0.0013
    ghbar = 0 #5e-6 
    catau = 50
    gcanbar = 0 #1e-7
    
    # Channel distribution
    Nadendfac = 0.035
    Kdendfac =  0.015
    ihdendfac = 1
    ldendfac = 0.25
    iahpdendfac = 0.1
    itinc = 2.39/60
    icaninc = itinc

    ###################################################################
    ## INSERT STIMULATION ELECTRODES
    stim = h.IClamp(.5)
    stim.delay = 1000
    stim.dur = dur
    stim.amp = istim

    ###################################################################
    ## INSERT ION CHANNELS:
    for sec in h.allsec():
        sec.insert("pas")
        sec.e_pas = Vrest
        sec.g_pas = 1/Rm
        sec.Ra = rall
        sec.cm = cap
        sec.insert("iar")
        sec.ghbar_iar = ghbar * ihdendfac
        sec.insert("Cad")
        sec.insert("ical")
        sec.insert("it2")
        sec.insert("iahp")
        sec.insert("hh2")
        sec.ena = 50
        sec.ek = -90
        sec.insert("ican")
     
     ##################################################################
    ## INSERT SYNAPSES    
    syn = h.Exp2Syn(h.dend[den](1)) # Sum of exp's
    syn.tau1 = 0.5 # Rise
    syn.tau2 = 2 # Decay
    syn.e = 10 # Reversal # pot.
              
    s = h.NetStim(0.5)
    s.start = 1000  # start for distal synapses
    s.number = 1
    s.noise = 0    
    
    nc = h.NetCon(s,syn,sec = sec)       
    nc.weight[0] = dsynweight

    #################################################################
    # for sec in h.soma:
    freq = 50
    h.geom_nseg(freq)
    tot = 0
    for sec in h.allsec():
        tot += sec.nseg
    h.distance()
    print("total # of segments (50Hz):", tot)      
              
    ##################################################################
    # INITIALIZE
    def initialize(Tchoice):
        global Epas
        h.celsius = celsius
        for sec in h.soma:
            h.distance()    
   
        for sec in h.allsec():
            sec.v = Vrest
            sec.e_pas = Epas
    
            sec.insert("pas")
            sec.e_pas = Epas
            sec.g_pas = 1/Rm
            sec.Ra = rall
            sec.cm = cap

        for sec in h.allsec():
            sec.gnabar_hh2 = gna * Nadendfac
            sec.vtraubNa_hh2 = nash  
            sec.gkbar_hh2 = gkdr * Kdendfac
            sec.vtraubK_hh2 = kdrsh
            sec.pcabar_ical = gcal * ldendfac
            sec.gkbar_iahp = gahp * iahpdendfac
            sec.ghbar_iar = ghbar * ihdendfac

        for sec in h.soma:
            sec.gnabar_hh2 = gna 
            sec.vtraubNa_hh2 = nash
            sec.gkbar_hh2 = gkdr
            sec.vtraubK_hh2 = kdrsh
            sec.pcabar_ical = gcal
            sec.gkbar_iahp = gahp
            sec.ghbar_iar = ghbar       
# Insert CaT
            sec.gcabar_it2 = gcat
            sec.gbar_ican = gcanbar            
            
        if (Tchoice == 7):  # Halnes
            Epas = -70.61
            for sec in h.allsec():
                sec.gcabar_it2 = gcat * (1 + itinc * h.distance(1))
                sec.gbar_ican = gcanbar * (1 + itinc * h.distance(1))
                
            for sec in h.soma:
                sec.gcabar_it2 = gcat
                sec.gbar_ican = gcanbar   
                  
        elif (Tchoice == 5): # Linear
              Epas = -70.6
              for sec in h.allsec():
                  for seg in sec:
                      seg.gcabar_it2 = gcat * (1 + itinc * (h.distance(0) + sec.L * seg.x)) 
                      seg.gbar_ican = gcanbar * (1 + itinc * (h.distance(0) + sec.L * seg.x))
      
        elif (Tchoice == 3): # Uniform
              Epas = -70.66
              for sec in h.allsec():        
                  sec.gcabar_it2 = gcat
                  sec.gbar_ican = gcanbar
                                              
        elif (Tchoice == 2): # Proximal
              Epas = -70.77
              mu = 60 # 
              sigma = 41 # Chosen to match data (Munsch et al, 1997) for 1st 60 micrometres.

              for sec in h.allsec():
                  for seg in sec:
                      d = h.distance(0) + sec.L * seg.x 
                      seg.gcabar_it2 = gcat * norm.pdf(d,mu,sigma)
                      seg.gbar_ican = gcanbar * norm.pdf(d,mu,sigma)
                      
        elif (Tchoice == 6): # Distal
              Epas = -70.72
              mu = 660 # 
              sigma = 120 #

              for sec in h.allsec():
                  for seg in sec:
                      d = h.distance(0) + sec.L * seg.x 
                      seg.gcabar_it2 = gcat * norm.pdf(d,mu,sigma)
                      seg.gbar_ican = gcanbar * norm.pdf(d,mu,sigma)              
                      
        elif (Tchoice == 4): # Middle
              Epas = -70.64
              mu = 300 # 
              sigma = 80 #

              for sec in h.allsec():
                  for seg in sec:
                      d = h.distance(0) + sec.L * seg.x 
                      seg.gcabar_it2 = gcat * norm.pdf(d,mu,sigma)
                      seg.gbar_ican = gcanbar * norm.pdf(d,mu,sigma) 
        
        elif (Tchoice == 0): # Null
              Epas = -69 
              for sec in h.allsec():
                  sec.gcabar_it2 = 0
                  sec.gbar_ican = 0                                                         
                                                                               
        elif (Tchoice == 1): # Soma
              Epas = -70.82   
              for sec in h.dend:
                  sec.gcabar_it2 = 0
                  sec.gbar_ican = 0    
        
        for sec in h.allsec():
            sec.taur_Cad = catau
           
        h.finitialize(Vrest)       
        h.fcurrent()     
        cvode.re_init()
    
    ###################################################################
    # Separate function that determines gcat
    def tcount():
        GCATOT = 0
        for sec in h.allsec():
            for seg in sec:    
                GCATOT += seg.gcabar_it2 * h.area(seg.x)  
        return GCATOT    
    
    ##################################################################
    # NORMALIZE g_CaT
    # Aim: Have the tot. # of T-channels as in Halnes et al. 2011 regardless dist. chosen.  
    initialize(7)    
    GCATOT7 = tcount()   
    print('GCATOT7', GCATOT7)
    initialize(Tdist)
    GCATOTX = tcount()
    print('GCATOTX', GCATOTX)
    if (Tdist !=0):
        gcat *=  GCATOT7/GCATOTX

    initialize(Tdist)

    print("gcat =", gcat)
    print("stim.amp=",stim.amp)

    ##################################################################
    vec ={}
    for var in 't', 'd_sec', 'd_seg', 'diam_sec','gc','diam_seg','stim_curr':
        vec[var] = h.Vector()
        
    for var in 'V_sec', 'V_seg', 'CaConc_sec','CaConc_seg':
        vec[var] = h.List()
    
    def create_lists(vec):        
        for sec in h.allsec():
            vec['d_sec'].append(h.distance(1))
            vec['diam_sec'].append(sec.diam)  
            rec0 = h.Vector()
            rec0.record(sec(0.5)._ref_v)         
            vec['V_sec'].append(rec0)
            rec_Ca = h.Vector()
            rec_Ca.record(sec(0.5)._ref_Cai)
            vec['CaConc_sec'].append(rec_Ca)        
            for seg in sec: 
                vec['d_seg'].append(h.distance(0) + sec.L * seg.x)
                vec['diam_seg'].append(seg.diam)
                vec['gc'].append(seg.gcabar_it2)   
                rec = h.Vector()
                rec.record(seg._ref_v)
                vec['V_seg'].append(rec)
                rec1 = h.Vector()
                rec1.record(seg._ref_Cai)
                vec['CaConc_seg'].append(rec1)                                   
        return vec
    
    #####################################    	
    create_lists(vec)       
    # run the simulation
    
    vec['t'].record(h._ref_t)
    # vec['current'].record(VC_patch._ref_i)
    vec['stim_curr'].record(stim._ref_i)
    h.load_file("stdrun.hoc")               
    h.tstop = 2500  # Simulation time
    h.t = -500
    h.run()
    return vec            

vec = test(Tdist,dur,gna,istim,dsynweight)

print("Dist #", Tdist,":",dist[Tdist], "; Epas =", Epas)

t = vec['t'].to_python()
d_seg = vec['d_seg'].to_python()
d_sec = vec['d_sec'].to_python()
gc = vec['gc'].to_python()
V_soma = vec['V_sec'][0].to_python()
CaConc_soma = vec['CaConc_sec'][0].to_python()
stim_curr = vec['stim_curr'].to_python()

# ########################################################################
## Plotting propagation of voltage signal 
fig = pylab.figure()

ax1 = fig.add_subplot(4, 1, 1)
ax1.set_title(dist[Tdist])
ax1.text(-0.025, 1.025, 'A',horizontalalignment='center',verticalalignment='bottom',fontsize=18, fontweight='demibold',
transform=ax1.transAxes)

ax2 = fig.add_subplot(4, 1, 2, sharex=ax1, sharey=ax1, ylabel='Voltage(mV)')
ax2.text(-0.025, 1.025, 'B',horizontalalignment='center',verticalalignment='bottom',fontsize=18, fontweight='demibold',
transform=ax2.transAxes)

ax3 = fig.add_subplot(4, 1, 3, sharex=ax1, sharey=ax1)
ax3.text(-0.025, 1.025, 'C',horizontalalignment='center',verticalalignment='bottom',fontsize=18, fontweight='demibold',
transform=ax3.transAxes)

ax4 = fig.add_subplot(4, 1, 4, sharex=ax1, sharey=ax1)
ax4.text(-0.025, 1.025, 'D',horizontalalignment='center',verticalalignment='bottom',fontsize=18, fontweight='demibold',
transform=ax4.transAxes)

ax1.plot(vec['t'], vec['V_sec'][0],color=colour[dist[Tdist]], label = 'At soma') 
ax2.plot(vec['t'], vec['V_sec'][83],color=colour[dist[Tdist]], label = '%1.1f $\mu$m' % d_sec[83]) 
ax3.plot(vec['t'], vec['V_sec'][88],color=colour[dist[Tdist]], label = '%1.1f $\mu$m' % d_sec[88])
ax4.plot(vec['t'], vec['V_sec'][99],color=colour[dist[Tdist]], label = '%1.1f $\mu$m' % d_sec[99]) 
try:
    ax1.legend(loc='best',frameon=False, fontsize=13)
    ax2.legend(loc='best',frameon=False, fontsize=13)
    ax3.legend(loc='best',frameon=False, fontsize=13)
    ax4.legend(loc='best',frameon=False, fontsize=13)
except:
    # no fontsize in older version of matplotlib
    ax1.legend(loc='best',frameon=False)
    ax2.legend(loc='best',frameon=False)
    ax3.legend(loc='best',frameon=False)
    ax4.legend(loc='best',frameon=False)

pylab.setp([a.set_xlabel('') for a in [ax1,ax2,ax3]], visible=False) 
pylab.setp([a.get_xticklabels() for a in [ax1,ax2,ax3]], visible=False)

    
