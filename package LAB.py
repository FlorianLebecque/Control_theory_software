def PID_RT(SP, PV, Man, MVMan, MVFF, Kc, Ti, Td, alpha, Ts, MVMin, MVMax, MV, MVP, MVI, MVD, E, ManFF=False, PVInit = 0, method = 'EBD-EBD') : 
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
    MVP.append(E[-1] * Kc)
    
    
    



    '''Mode manuel et anti wind-up'''
    if ManFF:
        MVFFI = MVFF[-1]
    else:
        MVFFI = 0
        
    if bool(Man[-1]) == True:
        if ManFF:
            MVI[-1] = MVMan[-1] - MVP[-1] - MVD[-1]
        else:
            MVI[-1] = MVMan[-1] - MVP[-1] - MVD[-1] - MVFFI

    '''Limitation de MV'''

    MV_TEMP = MVP[-1] + MVI[-1] + MVD[-1] + MVFFI
    
    if MV_TEMP >= MVMax:
        MVI[-1] = MVMax - MVP[-1] - MVD[-1] - MVFFI
        MV_TEMP = MVMax
        
    if MV_TEMP <= MVMin:
        MVI[-1] = MVMin - MVP[-1] - MVD[-1] - MVFFI
        MV_TEMP = MVMin
                   
    MV.append(MV_TEMP)
        

def LeadLag_RT(MV,Kp,T_lead,T_lag,Ts,PV,PVInit=0,method='EBD'):
    
    """
    The function "LeadLag_RT" needs to be included in a "for or while loop".
    
    :MV: input vector
    :Kp: process gain
    :T: lag time constant [s]
    :Ts: sampling period [s]
    :PV: output vector
    :PVInit: (optional: default value is 0)
    :method: discretisation method (optional: default value is 'EBD')
        EBD: Euler Backward difference
        EFD: Euler Forward difference
        TRAP: Trapezoïdal method
    
    The function "LeadLag_RT" appends a value to the output vector "PV".
    The appended value is obtained from a recurrent equation that depends on the discretisation method.
    """    
    if (T_lag != 0):
        K = Ts/T_lag
        if len(PV) == 0:
            PV.append(PVInit)
        else: # MV[k+1] is MV[-1] and MV[k] is MV[-2]
            if method == 'EBD':
                PV.append((1/(1+K)) * PV[-1] + ((K*Kp)/(1+K)) * ((1+(T_lead/Ts)) * MV[-1] - (T_lead/Ts) * MV[-2]))
            elif method == 'EFD':
                PV.append((1-K) * PV[-1] + K * Kp * ((T_lead/Ts) * MV[-1] + (1 - (T_lead/Ts)) * MV[-2]))
            elif method == 'TRAP':
                PV.append((1/(2 * T_lag+Ts)) * ((2 * T_lag-Ts) * PV[-1]+(2 * T_lead+Ts) * Kp * MV[-1] + (Ts-2*T_lead) * Kp * MV[-2]))
            else:
                PV.append((1/(1+K)) * PV[-1] + ((K*Kp)/(1+K)) * ((1+(T_lead/Ts)) * MV[-1] - (T_lead/Ts) * MV[-2]))
    else:
        PV.append(Kp*MV[-1])