def PID_RT(SP, PV, MAN, MVMan, MVFF, Kc, Ti, Td, alpha, Ts, MVMin, MVMax, MV, MVP, MVI, MVD, E, ManFF=False, PVInit = 0, method = 'EBD-EBD') : 
    '''
The function "PID_RT" needs to be included in a "for or while loop". 

:SP: SP (or SetPoint) vector 
:PV: PV (or Process Value) vector 
:Man: Man (or Manual controller mode) vector (True or False) 
:MVMan: MVMan (or Manual value for MV) vector 
:MVFF: MVFF (or Feedforward) vector 

:Kc: controller gain 
:Ti: integral time constant [s] 
:Td: derivative time constant [s] 
:alpha: Tfd = alpha*Td where Tfd is the derivative filter time constant [s] 
:Ts: sampling period [s] 

:MVMin: minimum value for MV (used for saturation and anti wind-up) 
:MVMax: maximum value for MV (used for saturation and anti wind-up) 

:MV: MV (or Manipulated Value) vector 
:MVP: MVP (or Propotional part of MV) vector 
:MVI: MVI (or Integral part of MV) vector 
:MVD: MVD (or Derivative part of MV) vector 
:E: E (or control Error) vector 

:ManFF: Activated FF in manual mode (optional: default boolean value is False) 
:PVInit: Initial value for PV (optional: default value is 0): used if PID_RT is ran first in the squence and no value of PV is available yet. 

:method: discretisation method (optional: default value is 'EBD') 
    EBD-EBD: EBD for integral action and EBD for derivative action 
    EBD-TRAP: EBD for integral action and TRAP for derivative action 
    TRAP-EBD: TRAP for integral action and EBD for derivative action 
    TRAP-TRAP: TRAP for integral action and TRAP for derivative action 

The function "PID_RT" appends new values to the vectors "MV", "MVP", "MVI", and "MVD". The appended values are based on the PID algorithm, the controller mode, and feedforward. Note that saturation of "MV" within the limits [MVMin MVMax] is implemented with anti wind-up. 
    '''
    if len(PV) == 0 : 
        E.append(SP[-1] - PVInit)
    else : 
        E.append(SP[-1] - PV[-1])
    
    methodI, methodD =  method.split('-')
    
    
    if len(MVI) == 0 : 
        MVI.append((Kc*Ts/Ti)*E[-1])
    else : 
        if methodI == 'TRAP':
            MVI.append(MVI[-1] + (0.5*Kc*Ts/Ti)*(E[-1]+E))
        else : 
            MVI.append(MVI[-1] + (Kc*Ts/Ti)*E[-1])
            
    '''Initialisation de MVD : à 0 directement, voir slide 193  '''
    Tfd = alpha*Td
    if len(MVD) == 0 : 
        MVD.append(0)
    else : 
        if methodD == 'EBD':
            MVD.append( ((Tfd)/(Tfd + Ts))*MVD[-1] + ((Kc*Td)/(Tfd+Ts))*E[-1] - E[-2])
        
    '''Initialisation de MVP'''
    MVP.append(E[-1]*Kc)
    
    
    ''' Application de MVP, MVK et MVI sur MV, slide 194'''
    MV.append(MVP[-1] + MVK[-1] + MVI[-1])
    
        