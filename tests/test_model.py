# model1: DDC + MAO pathway, no DHPAA drain, +product inhibition
def model1(u, t, k0, Vmax1, KmA, Vmax2, KmB, Kmo2, O2, k3, k4, Inh, KiB, KiC):
    
    A = u[0] # Dopa
    B = u[1] # Dopamin
    C = u[2] # DHPAA
    D = u[3] # THP
    S = u[4] # Total fed Dopa
    P = u[5] # Total converted THP
    
    dAdt = k0 * (1.0 - pow(S/100, 30)) - Vmax1 * A / (KmA * (1 + B / KiB) + A)
    dBdt = Vmax1 * A / (KmA * (1 + B / KiB) + A) - Vmax2 * B * O2 / ((KmB * (1 + Inh + C / KiC) + B) * (Kmo2 + O2)) - k3 * B * C
    dCdt = Vmax2 * B * O2 / ((KmB * (1 + Inh + C / KiC) + B) * (Kmo2 + O2)) - k3 * B * C
    dDdt = k3 * B * C - k4 * D
    dSdt = k0 * (1.0 - pow(S/100, 30))
    dPdt = k3 * B * C
    return([dAdt, dBdt, dCdt, dDdt, dSdt, dPdt])
    
    
    
# model2: DDC + MAO pathway, +DHPAA drain, +product inhibition
def model2(u, t, k0, Vmax1, KmA, Vmax2, KmB, Kmo2, O2, k3, k4, Inh, KiB, KiC, k5):
    
    A = u[0] # Dopa
    B = u[1] # Dopamin
    C = u[2] # DHPAA
    D = u[3] # THP
    S = u[4] # Total fed Dopa
    P = u[5] # Total converted THP
    
    dAdt = k0 * (1.0 - pow(S/100, 30)) - Vmax1 * A / (KmA * (1 + B / KiB) + A)
    dBdt = Vmax1 * A / (KmA * (1 + B / KiB) + A) - Vmax2 * B * O2 / ((KmB * (1 + Inh + C / KiC) + B) * (Kmo2 + O2)) - k3 * B * C
    dCdt = Vmax2 * B * O2 / ((KmB * (1 + Inh + C / KiC) + B) * (Kmo2 + O2)) - k3 * B * C - k5 * C
    dDdt = k3 * B * C - k4 * D
    dSdt = k0 * (1.0 - pow(S/100, 30))
    dPdt = k3 * B * C
    return([dAdt, dBdt, dCdt, dDdt, dSdt, dPdt])

# model3: DDC + DHPAAS pathway, no DHPAA drain, +product inhibition
def model3(u, t, k0, Vmax1, KmA, Vmax2, KmA2, Kmo2, O2, k3, k4, KiB, KiC):
    
    A = u[0] # Dopa
    B = u[1] # Dopamin
    C = u[2] # DHPAA
    D = u[3] # THP
    S = u[4] # Total fed Dopa
    P = u[5] # Total converted THP
    
    dAdt = k0 * (1.0 - pow(S/100, 30)) - Vmax1 * A / ((KmA * (1 + B / KiB) + A)) - Vmax2 * A * O2 / ((KmA2 * (1 + C / KiC) + A) * (Kmo2 + O2))
    dBdt = Vmax1 * A / ((KmA * (1 + B / KiB) + A)) - k3 * B * C
    dCdt = Vmax2 * A * O2 / ((KmA2 * (1 + C / KiC) + A) * (Kmo2 + O2)) - k3 * B * C
    dDdt = k3 * B * C - k4 * D
    dSdt = k0 * (1.0 - pow(S/100, 30))
    dPdt = k3 * B * C
    return([dAdt, dBdt, dCdt, dDdt, dSdt, dPdt])

# model4: DDC + DHPAAS pathway, +DHPAA drain, +product inhibition
def model4(u, t, k0, Vmax1, KmA, Vmax2, KmA2, Kmo2, O2, k3, k4, KiB, KiC, k5):
    
    A = u[0] # Dopa
    B = u[1] # Dopamin
    C = u[2] # DHPAA
    D = u[3] # THP
    S = u[4] # Total fed Dopa
    P = u[5] # Total converted THP
    
    dAdt = k0 * (1.0 - pow(S/100, 30)) - Vmax1 * A / ((KmA * (1 + B / KiB) + A)) - Vmax2 * A * O2 / ((KmA2 * (1 + C / KiC) + A) * (Kmo2 + O2))
    dBdt = Vmax1 * A / ((KmA * (1 + B / KiB) + A)) - k3 * B * C
    dCdt = Vmax2 * A * O2 / ((KmA2 * (1 + C / KiC) + A) * (Kmo2 + O2)) - k3 * B * C - k5 * C
    dDdt = k3 * B * C - k4 * D
    dSdt = k0 * (1.0 - pow(S/100, 30))
    dPdt = k3 * B * C
    return([dAdt, dBdt, dCdt, dDdt, dSdt, dPdt])

