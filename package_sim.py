


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
                value = (1 / 1 + K) * PV[-1] + (((Kp * K)/(1 + K)) * (((1 + (T_lead / Ts)) * MV[-1]) - ((T_lead / Ts) * MV[-2])))
            elif method == 'EFD':
                value = ((1-K) * PV[-1]) + ((Kp * K) * ( (T_lead/Ts) * MV[-1] ) + ((1 - (T_lead / Ts)) * MV[-2]) )
            else:
                value = (1 / 1 + K) * PV[-1] + (((Kp * K)/(1 + K)) * (((1 + (T_lead / Ts)) * MV[-1]) - ((T_lead / Ts) * MV[-2])))

            PV.append(value)
            
    else:
        PV.append(Kp*MV[-1])