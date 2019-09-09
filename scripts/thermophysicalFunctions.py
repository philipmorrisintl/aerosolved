import numpy as np

def psPG(T):
    return max(np.exp(212.8 - 15420/T - 28.109 * np.log(T) + 2.1564E-5*T**2), 0)

def psVG(T):
    return 7.5E6 * np.exp \
    ( \
        (850.0/min(T,850)) \
      * ( \
          - 6.94758*max(1.0 - min(T,850.0)/850.0, 0.0) \
          - 0.33345*max(1.0 - min(T,850.0)/850.0, 0.0)**(3.0/2.0) \
          - 5.98569*max(1.0 - min(T,850.0)/850.0, 0.0)**(5.0/2.0) \
          - 1.33011*max(1.0 - min(T,850.0)/850.0, 0.0)**5 \
        ) \
    )

def psWater(T):
    max(np.exp(73.649 - 7258.2/T - 7.3037 * np.log(T) + 4.1653E-6*T**2), 0)

def rhoPGl(T):
    return max \
    ( \
        83.11748 \
      / np.power \
        ( \
            0.26106, \
            1.0 \
          + np.power \
            ( \
               min(max(1.0- min(T/626,1.0),1E-16), 1.0), \
               0.20459 \
            ) \
        ), \
        0 \
    )

def rhoVGl(T):
    return 349.0 \
  + 1341.5932 * min(max(1.0-T/850.0, 0.0), 1.0)**0.35 \
  - 1168.205  * min(max(1.0-T/850.0, 0.0), 1.0)**(2.0/3.0) \
  + 1429.7634 * min(max(1.0-T/850.0, 0.0), 1.0) \
  - 527.771   * min(max(1.0-T/850.0, 0.0), 1.0)**(4.0/3.0)

def rhoWaterl(T):
    return max(((3.280712e-05*T - 0.03440865)*T + 11.53645)*T -249.5258, 0)

def rhoPGg(T):
    return 76.094/8.3144621/T*1E2

def rhoVGg(T):
    return 92.09/8.3144621/T*1E2

def rhoWaterg(T):
    return 18.015/8.3144621/T*1E2

def rhoAirg(T):
    return 28.810/8.3144621/T*1E2

def muPGg(T):
    As = 1.67212e-06
    Ts = 170.672
    return As*np.sqrt(T)/(1.0+Ts/T)

def muVGg(T):
    As = 1.67212e-06
    Ts = 170.672
    return As*np.sqrt(T)/(1.0+Ts/T)

def muWaterg(T):
    As = 1.67212e-06
    Ts = 170.672;
    return As*np.sqrt(T)/(1.0+Ts/T)

def muAirg(T):
    As = 1.67212e-06
    Ts = 170.672;
    return As*np.sqrt(T)/(1.0+Ts/T)

