


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

                TLeadDivTs = T_lead / Ts
                OneDivOnePlusK = (1 / (1+K))
                mult_part = ((Kp*K)/(1+K))
                one_plus_tleadDivTs = 1 + TLeadDivTs

                value = (OneDivOnePlusK * PV[-1]) + ( mult_part * ( (one_plus_tleadDivTs * MV[-1]) - (TLeadDivTs * MV[-2] )) )

            elif method == 'EFD':
                OneMinusK = 1-K
                KTimesKp = K*Kp
                TLeadDivTs = T_lead / Ts
                one_minus_tleadDivTs = 1 - TLeadDivTs

                value = (OneMinusK * PV[-1]) + ( KTimesKp * ( (TLeadDivTs * MV[-1]) + (one_minus_tleadDivTs * MV[-2]) ) )
            else:
                TLeadDivTs = T_lead / Ts
                OneDivOnePlusK = (1 / (1+K))
                mult_part = ((Kp*K)/(1+K))
                one_plus_tleadDivTs = 1 + TLeadDivTs

                value = (OneDivOnePlusK * PV[-1]) + ( mult_part * ( (one_plus_tleadDivTs * MV[-1]) - (TLeadDivTs * MV[-2] )) )

            PV.append(value)
            
    else:
        PV.append(Kp*MV[-1])