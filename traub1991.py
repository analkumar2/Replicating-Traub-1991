#Author: Anal Kumar
#Uses the parameters in Traub et al, 1991 to model a reduced CA3 neuron
#Morphology was already set up by Harshavardan


import moose
import rdesigneur as rd
import numpy as np
import matplotlib.pyplot as plt

EREST = -0.060
F = 96485.3329

#R is the radius of the compartment. B is the scaling term.
#surfA is the surface area of the compartment.
#The function finds the thickness of the Ca shell when the scaling
#term is given.
def findThickness(R, surfA, B):
    thick = R - np.sqrt(R**2 - 2*R/(B*F*surfA))
    return(str(thick))



rdes = rd.rdesigneur(
    cellProto = [
        ['Compartments.swc','elec']
    ],
        
    chanProto = [
        ['Channelprotos.NaChan()', 'Na_Chan'],
        ['Channelprotos.KdrChan()', 'Kdr_Chan'],
        ['Channelprotos.KaChan()', 'Ka_Chan'],
        ['Channelprotos.CaConc()', 'Ca_conc'],
        ['Channelprotos.CaChan()', 'Ca_Chan'],        
        ['Channelprotos.KahpChan()', 'Kahp_Chan'],
        ['Channelprotos.KcChan()', 'Kc_Chan'],
        ['Channelprotos.LChan()', 'L_Chan'],
    ],
        
    passiveDistrib = [
        ['soma', 'Rm', '301204819.3', 'Ra', '37650.60241', 'Cm', '9.96E-11', 'initVm', str(EREST), 'Em', str(EREST)],
        ['apical#', 'RM', '1', 'RA', '1', 'CM', '0.03', 'initVm', str(EREST), 'Em', str(EREST)],
        ['dend#', 'RM', '1', 'RA', '1', 'CM', '0.03', 'initVm', str(EREST), 'Em', str(EREST)],
    ],
        
    chanDistrib = [
        #soma
        ['Na_Chan', 'soma', 'Gbar', '300'],
        ['Kdr_Chan', 'soma', 'Gbar', '150'],
        ['Ka_Chan', 'soma', 'Gbar', '50'],
        ['Ca_conc', 'soma', 'thick', findThickness(16.25e-6, 3320e-12, 17.402e12)],
        ['Ca_Chan', 'soma', 'Gbar', '40'],        
        ['Kahp_Chan', 'soma', 'Gbar', '8'],
        ['Kc_Chan', 'soma', 'Gbar', '100'],
        ['L_Chan', 'soma', 'Gbar', '1'],
        
        # apical dendrtites
        ['Na_Chan', 'apical#', 'Gbar', '(p<=120e-6) ? 150 : ((p>=240e-6 && p<=360e-6) ? 200 : 0)'],
        ['Kdr_Chan', 'apical#', 'Gbar', '(p<=120e-6) ? 50 : ((p>=240e-6 && p<=360e-6) ? 200 : 0)'],
        ['Ca_conc', 'apical#', 'thick', '(p<=120e-6) ? %s : %s' %(findThickness(2.89e-6, 2188e-12, 26.404e12), findThickness(2.89e-6, 2188e-12, 5.941e12))],
        ['Ca_Chan', 'apical#', 'Gbar', '(p<=120e-6) ? 80 : ((p>=120e-6 && p<=240e-6) ? 50 : ((p>=240e-6 && p<=600e-6) ? 170 : ((p>=600e-6 && p<=840e-6) ? 100 : ((p>=840e-6 && p<=1080e-6) ? 50 : 0))))'],
        ['Kahp_Chan', 'apical#', 'Gbar', '(p<=1080e-6) ? 8 : 0'],
        ['Kc_Chan', 'apical#', 'Gbar', '(p<=120e-6) ? 200 : ((p>=120e-6 && p<=240e-6) ? 50 : ((p>=240e-6 && p<=840e-6) ? 100 : ((p>=840e-6 && p<=1080e-6) ? 50 : 0)))'],
        ['L_Chan', 'apical#', 'Gbar', '1'],
        
        # basal dendrites
        ['Na_Chan', 'dend#', 'Gbar', '(p<=110e-6) ? 150 : ((p>=220e-6 && p<=330e-6) ? 200 : 0)'],
        ['Kdr_Chan', 'dend#', 'Gbar', '(p<=110e-6) ? 50 : ((p>=220e-6 && p<=330e-6) ? 200 : 0)'],
        ['Ca_conc', 'dend#', 'thick', '(p<=110e-6) ? %s : %s' %(findThickness(2.42e-6, 1673e-12, 34.530e12), findThickness(2.42e-6, 1673e-12, 7.769e12))],
        ['Ca_Chan', 'dend#', 'Gbar', '(p<=110e-6) ? 80 : ((p>=110e-6 && p<=220e-6) ? 50 : ((p>=220e-6 && p<=550e-6) ? 120 : ((p>=550e-6 && p<=770e-6) ? 50 : 0)))'],
        ['Kahp_Chan', 'dend#', 'Gbar', '(p<=770e-6) ? 8 : 0'],
        ['Kc_Chan', 'dend#', 'Gbar', '(p<=110e-6) ? 200 : ((p>=110e-6 && p<=220e-6) ? 50 : ((p>=220e-6 && p<=550e-6) ? 100 : ((p>=550e-6 && p<=770e-6) ? 50 : 0)))'],
        ['L_Chan', 'dend#', 'Gbar', '1'],
    ],
        
    stimList = [
        ['soma', '1', '.', 'inject', '(t>1 && t<=1.2) ? 1e-8 : 0' ],
    ],
    
    plotList = [
        ['soma', '1', '.', 'Vm', 'Membrane potential'],
        ['soma', '1', 'Ca_Chan', 'Ik', 'Calcium current'],
        ['soma', '1', 'Ca_conc', 'Ca', 'Calcium concentration'],
    ]
        
    # moogList = [
        # ['#', '1', '.', 'Vm', 'Soma potential'],
    # ]
)

rdes.buildModel()
moose.reinit()
moose.start( 2 )
rdes.display()