 PARAM_DIST = {
                "EGP_0": {"mean": 0.0169, "std": 0.0039, "unit": "mmol/min", "distribution": "Gaussian", "variability": [[0, 3], 0.05]},
                "F_01": {"mean": 0.0111, "std": 0.0007, "unit": "mmol/min", "distribution": "Gaussian", "variability": [[0, 3], 0.05]},
                "k_12": {"mean": 0.00649, "std": 0.00282, "unit": "1/min", "distribution": "Gaussian", "variability": [[0, 3], 0.05]},
                "S_id": {"mean": 0.00082, "std": 0.00032, "unit": "mU/(L min)", "distribution": "Gaussian", "variability": [[0, 3], 0.05]},
                "S_ie": {"mean": 0.052, "std": 0.0125, "unit": "L/mU", "distribution": "Gaussian", "variability": [[0, 3], 0.05]},
                "S_it": {"mean": 0.00512, "std": 0.00131, "unit": "mU/(L min)", "distribution": "Gaussian", "variability": [[0, 3], 0.05]},
                "k_e": {"mean": 0.14, "std": 0.035, "unit": "1/min", "distribution": "Gaussian", "variability": None},
                "Bio": {"mean": 0.7, "std": 1.2, "unit": "%", "distribution": "Uniform", "variability": [[0, 24], 0.2]},
                "k_a_int": {"mean": -2.372, "std": 1.092, "unit": "1/min", "distribution": "Lognormal", "variability": [[0, 3], 0.05]},
                "k_a1": {"mean": 0.0055, "std": 0.0056, "unit": "1/min", "distribution": "Gaussian", "variability": [[0, 3], 0.05]},
                "k_a2": {"mean": 0.0683, "std": 0.0507, "unit": "1/min", "distribution": "Gaussian", "variability": [[0, 3], 0.05]},
                "k_a3": {"mean": 0.0304, "std": 0.0235, "unit": "1/min", "distribution": "Gaussian", "variability": [[0, 3], 0.05]},
                "V_g": {"mean": -1.8971199848858813, "std": 0.23, "unit": "L/kg", "distribution": "Lognormal", "variability": None},
                "V_i": {"mean": 0.12, "std": 0.012, "unit": "L/kg", "distribution": "Gaussian", "variability": None},
                "k_a": {"mean": 0.018, "std": 0.0045, "unit": "mmol/L", "distribution": "Gaussian", "variability": [[0, 3], 0.05]},
                "t_max": {"mean": -3.689, "std": 0.025, "unit": "min", "distribution": "Lognormal", "variability": None},
                "w": {"mean": 74.9, "std": 14.4, "unit": "kg", "distribution": "Gaussian", "variability": None},
                "R_th":{"mean": 9, "std": 1.5, "unit": "mmol/L", "distribution": "Gaussian", "variability": None},
                "R_cl": {"mean": 0.01, "std": 0.025, "unit": "1/min", "distribution": "Gaussian", "variability": None},
                "U_ceil": {"mean": 0.02, "std": 0.035, "unit": "TODO ??", "distribution": "Uniform", "variability": None}
        }

old_data
EGP_0_dist = [0.0169, 0.0039]   # mmol / min
        F_01_dist = [0.0111, 0.0007]    # mmol / min
        k_12_dist = [0.00649, 0.00282]  # 1 / min 
        S_it_dist = [0.00512, 0.00131]  # 1 / min (S_I1)
        S_id_dist = [0.00082, 0.00032]  # 1 / min (S_I2)
        S_ie_dist = [0.052, 0.0125]     # 1 / min (S_I3)
        k_e_dist = [0.14, 0.035]        # 1 / min
        Bio_dist = [.70, 1.20]          #
        k_aint_dist = [-2.372, 1.092]   # 1 / min (1 / T_I)

# Intra and Inter variation
        EGP_0_dist = [0.0169, 0.0039]   # mmol / min
        F_01_dist = [0.0111, 0.0007]    # mmol / min
        k_12_dist = [0.00649, 0.00282]  # 1 / min 
        S_it_dist = [0.00512, 0.00131]  # 1 / min (S_I1)
        S_id_dist = [0.00082, 0.00032]  # 1 / min (S_I2)
        S_ie_dist = [0.052, 0.0125]     # 1 / min (S_I3)
        k_e_dist = [0.14, 0.035]        # 1 / min
        Bio_dist = [.70, 1.20]          #
        k_aint_dist = [-2.372, 1.092]   # 1 / min (1 / T_I)

        # Intra individuals variation
        k_a1 = np.random.normal(0.0055, 0.0056)
        k_a2 = np.random.normal(0.0683, 0.0507)
        k_a3 = np.random.normal(0.0304, 0.0235)
        w = np.random.normal(74.9, 14.4)
        
        # Volumes depends on weight
        V_g = np.random.lognormal(np.log(0.15), 0.23) * w
        V_i = np.random.normal(0.12, .012) * w
        # V_g = np.log((np.random.normal(1.16, 0.23)))      
        
        k_a = np.random.normal(0.018, 0.0045)                # T_S
        t_max = 1 / np.exp(np.random.normal(-3.689, 0.025))  # T_D                     
        R_th = np.random.normal(9 ,1.5)
        R_cl = np.abs(np.random.normal(0.01 ,0.025))
        U_g_ceil = np.random.uniform(0.02, 0.035) * w

        # Oscillatory behavior of cyclic of variables: Sinusoidal Oscillation [U, z]   
        self.circadian_time = [[np.random.uniform(0,3), 0.05],  # EGP_0
                               [np.random.uniform(0,3), 0.05],  # F_01 
                               [np.random.uniform(0,3), 0.05],  # k_12
                               [np.random.uniform(0,3), 0.05],  # S_it 
                               [np.random.uniform(0,3), 0.05],  # S_id 
                               [np.random.uniform(0,3), 0.05],  # S_ie 
                               [np.random.uniform(0,3), 0.05],  # k_e 
                               [np.random.uniform(0,24), 0.2],  # Bio
                               [np.random.uniform(0,3), 0.05]]  # k_aint
        # Initial values of cyclic of variables [p0]   
        self.cyclic_params = np.array([ np.random.normal(EGP_0_dist[0], EGP_0_dist[1]), 
                                        np.random.normal(F_01_dist[0], F_01_dist[1]),
                                        np.random.normal(k_12_dist[0], k_12_dist[1]),
                                        np.random.normal(S_id_dist[0], S_id_dist[1]),
                                        np.random.normal(S_ie_dist[0], S_ie_dist[1]),
                                        np.random.normal(S_it_dist[0], S_it_dist[1]),
                                        np.random.normal(k_e_dist[0], k_e_dist[1]),
                                        np.random.uniform(Bio_dist[0], Bio_dist[1]),
                                        np.random.lognormal(k_aint_dist[0], k_aint_dist[1])], dtype='float64')
        
        # Parameters Initialization
        if params == None:
            # Inicialização dos Parâmetros Oscilatórios dos Pacientes
            # [   *0,   *1,   *2,   *3,   *4,   *5,  *6,  *7,     *8,    9,   10,  11,   12,  13,  14,    15,16,   17,   18,     19]
            # [EGP_0, F_01, k_12, S_id, S_ie, S_it, k_e, Bio, k_aint, k_a1, k_a2, k_a3, V_g, V_i, k_a, t_max, w, R_th, R_cl, U_ceil]




            ## ANTIGO #####

            # Intra and Inter variation
        EGP_0_dist = [0.0169, 0.0039]   # mmol / min
        F_01_dist = [0.0111, 0.0007]    # mmol / min
        k_12_dist = [0.00649, 0.00282]  # 1 / min 
        S_it_dist = [0.00512, 0.00131]  # 1 / min (S_I1)
        S_id_dist = [0.00082, 0.00032]  # 1 / min (S_I2)
        S_ie_dist = [0.052, 0.0125]     # 1 / min (S_I3)
        k_e_dist = [0.14, 0.035]        # 1 / min
        Bio_dist = [.70, 1.20]          #
        k_aint_dist = [-2.372, 1.092]   # 1 / min (1 / T_I)

        # Intra individuals variation
        k_a1 = np.random.normal(0.0055, 0.0056)
        k_a2 = np.random.normal(0.0683, 0.0507)
        k_a3 = np.random.normal(0.0304, 0.0235)
        w = np.random.normal(74.9, 14.4)
        
        # Volumes depends on weight
        V_g = np.random.lognormal(np.log(0.15), 0.23) * w
        V_i = np.random.normal(0.12, .012) * w
        # V_g = np.log((np.random.normal(1.16, 0.23)))      
        
        k_a = np.random.normal(0.018, 0.0045)                # T_S
        t_max = 1 / np.exp(np.random.normal(-3.689, 0.025))  # T_D                     
        R_th = np.random.normal(9 ,1.5)
        R_cl = np.abs(np.random.normal(0.01 ,0.025))
        U_g_ceil = np.random.uniform(0.02, 0.035) * w

        # Oscillatory behavior of cyclic of variables: Sinusoidal Oscillation [U, z]   
        self.circadian_time = [[np.random.uniform(0,3), 0.05],  # EGP_0
                               [np.random.uniform(0,3), 0.05],  # F_01 
                               [np.random.uniform(0,3), 0.05],  # k_12
                               [np.random.uniform(0,3), 0.05],  # S_it 
                               [np.random.uniform(0,3), 0.05],  # S_id 
                               [np.random.uniform(0,3), 0.05],  # S_ie 
                               [np.random.uniform(0,3), 0.05],  # k_e 
                               [np.random.uniform(0,24), 0.2],  # Bio
                               [np.random.uniform(0,3), 0.05]]  # k_aint
        # Initial values of cyclic of variables [p0]   
        self.cyclic_params = np.array([ np.random.normal(EGP_0_dist[0], EGP_0_dist[1]), 
                                        np.random.normal(F_01_dist[0], F_01_dist[1]),
                                        np.random.normal(k_12_dist[0], k_12_dist[1]),
                                        np.random.normal(S_id_dist[0], S_id_dist[1]),
                                        np.random.normal(S_ie_dist[0], S_ie_dist[1]),
                                        np.random.normal(S_it_dist[0], S_it_dist[1]),
                                        np.random.normal(k_e_dist[0], k_e_dist[1]),
                                        np.random.uniform(Bio_dist[0], Bio_dist[1]),
                                        np.random.lognormal(k_aint_dist[0], k_aint_dist[1])], dtype='float64')
        
        # Parameters Initialization
        if params == None:
            # Inicialização dos Parâmetros Oscilatórios dos Pacientes
            # [   *0,   *1,   *2,   *3,   *4,   *5,  *6,  *7,     *8,    9,   10,  11,   12,  13,  14,    15,16,   17,   18,     19]
            # [EGP_0, F_01, k_12, S_id, S_ie, S_it, k_e, Bio, k_aint, k_a1, k_a2, k_a3, V_g, V_i, k_a, t_max, w, R_th, R_cl, U_ceil]
            self.params = np.array([self.cyclic_params[0], 
                                    self.cyclic_params[1],
                                    self.cyclic_params[2],
                                    self.cyclic_params[3],
                                    self.cyclic_params[4],
                                    self.cyclic_params[5],
                                    self.cyclic_params[6],
                                    self.cyclic_params[7],
                                    self.cyclic_params[8],
                                    k_a1, k_a2, k_a3,
                                    V_g, V_i,
                                    k_a, t_max, w,
                                    R_th, R_cl, U_g_ceil], dtype='float64')   





















                                    
import matplotlib.pyplot as plt
import matplotlib.dates as md
from scipy.integrate import solve_ivp
from scipy.optimize import root
import numpy as np
import datetime

#################### Modelo Combinado de Oscilação ########################
#### TODO LIST
# - PARAMETERS
# - TRANSLATE
# - STEADY STATE u_s --OK?
# - CONTROL
# - 
class Subject():
    
    def __init__(self, params=None, food_intake_f=None, insulin_intake_f=None):
        """
        Class Subject
        
        > params(default None): list with all the parameters in Hovorka 's Model in order given bellow. 
        In case of None, automatic parameters will be given following the built-in probability distribution.
        [EGP_0, F_01, k_12, S_id, S_ie, S_it, k_e, Bio, k_a1, k_a2, k_a3, V_g, V_i, k_a, t_max, w, R_th, R_cl, U_ceil]
        
        > food_intake_f(default None): Function in format 'fun(t)' that return the rate of carbohydrate ingestetion in given time t.
        The values should be (g / min).
        
        > insulin_intake_f(default None): Function in format 'fun(t)' that return the rate of  insulin infusion in given time t.
        The values should be (U / min).

        """
        # Variable with resulting simulation values from 'simulate'
        self.solution = np.array([])
        self.time = np.array([])
        self.max_glucose = 0
        self.min_glucose = 0
        
        # Steady-state basal insulin value
        self.u_basal = 0

        # Auxiliary input functions: food intake and insulin infusion
        self.insulin_intake_fun = insulin_intake_f
        self.food_intake_fun = food_intake_f
        self.insulin = np.array([])        
        self.meal = np.array([])        
        
        # Sensor reading
        self.cgm_glucose = np.array([])
        self.cgm_time = np.array([])       

        # Paramete Distributions
        # Paper Software for in Silico Testing.pdf
        
        # Intra and Inter variation
        EGP_0_dist = [0.0169, 0.0039]   # mmol / min
        F_01_dist = [0.0111, 0.0007]    # mmol / min
        k_12_dist = [0.00649, 0.00282]  # 1 / min 
        S_it_dist = [0.00512, 0.00131]  # 1 / min (S_I1)
        S_id_dist = [0.00082, 0.00032]  # 1 / min (S_I2)
        S_ie_dist = [0.052, 0.0125]     # 1 / min (S_I3)
        k_e_dist = [0.14, 0.035]        # 1 / min
        Bio_dist = [.70, 1.20]          #
        k_aint_dist = [-2.372, 1.092]   # 1 / min (1 / T_I)

        # Intra individuals variation
        k_a1 = np.random.normal(0.0055, 0.0056)
        k_a2 = np.random.normal(0.0683, 0.0507)
        k_a3 = np.random.normal(0.0304, 0.0235)
        w = np.random.normal(74.9, 14.4)
        
        # Volumes depends on weight
        V_g = np.random.lognormal(np.log(0.15), 0.23) * w
        V_i = np.random.normal(0.12, .012) * w
        # V_g = np.log((np.random.normal(1.16, 0.23)))      
        
        k_a = np.random.normal(0.018, 0.0045)                # T_S
        t_max = 1 / np.exp(np.random.normal(-3.689, 0.025))  # T_D                     
        R_th = np.random.normal(9 ,1.5)
        R_cl = np.abs(np.random.normal(0.01 ,0.025))
        U_g_ceil = np.random.uniform(0.02, 0.035) * w

        # Oscillatory behavior of cyclic of variables: Sinusoidal Oscillation [U, z]   
        self.circadian_time = [[np.random.uniform(0,3), 0.05],  # EGP_0
                               [np.random.uniform(0,3), 0.05],  # F_01 
                               [np.random.uniform(0,3), 0.05],  # k_12
                               [np.random.uniform(0,3), 0.05],  # S_it 
                               [np.random.uniform(0,3), 0.05],  # S_id 
                               [np.random.uniform(0,3), 0.05],  # S_ie 
                               [np.random.uniform(0,3), 0.05],  # k_e 
                               [np.random.uniform(0,24), 0.2],  # Bio
                               [np.random.uniform(0,3), 0.05]]  # k_aint
        # Initial values of cyclic of variables [p0]   
        self.cyclic_params = np.array([ np.random.normal(EGP_0_dist[0], EGP_0_dist[1]), 
                                        np.random.normal(F_01_dist[0], F_01_dist[1]),
                                        np.random.normal(k_12_dist[0], k_12_dist[1]),
                                        np.random.normal(S_id_dist[0], S_id_dist[1]),
                                        np.random.normal(S_ie_dist[0], S_ie_dist[1]),
                                        np.random.normal(S_it_dist[0], S_it_dist[1]),
                                        np.random.normal(k_e_dist[0], k_e_dist[1]),
                                        np.random.uniform(Bio_dist[0], Bio_dist[1]),
                                        np.random.lognormal(k_aint_dist[0], k_aint_dist[1])], dtype='float64')
        
        # Parameters Initialization
        if params == None:
            # Inicialização dos Parâmetros Oscilatórios dos Pacientes
            # [   *0,   *1,   *2,   *3,   *4,   *5,  *6,  *7,     *8,    9,   10,  11,   12,  13,  14,    15,16,   17,   18,     19]
            # [EGP_0, F_01, k_12, S_id, S_ie, S_it, k_e, Bio, k_aint, k_a1, k_a2, k_a3, V_g, V_i, k_a, t_max, w, R_th, R_cl, U_ceil]
            self.params = np.array([self.cyclic_params[0], 
                                    self.cyclic_params[1],
                                    self.cyclic_params[2],
                                    self.cyclic_params[3],
                                    self.cyclic_params[4],
                                    self.cyclic_params[5],
                                    self.cyclic_params[6],
                                    self.cyclic_params[7],
                                    self.cyclic_params[8],
                                    k_a1, k_a2, k_a3,
                                    V_g, V_i,
                                    k_a, t_max, w,
                                    R_th, R_cl, U_g_ceil], dtype='float64')                      
        else:
            self.params = params
        
        # Steady-state basal insulin
        self.steady_state()
    
    def food_intake(self, t):
        """
        Food ingestion function. Default 50g at 30 minutes during 15 minutes.

        > t: time [min].
        < D: Equivalent glicose ingested [mmol /min].
        """

        if self.food_intake_fun == None:    
            # 50g during 15 minutes
            t0 = 30
            dt = 15
            d = 50 / dt if (t >= t0 and t <= t0 + dt) else 0
        else:
            d = self.food_intake_fun(t)

        # Convert from [g/min] to [mmol/min]
        return d * 1000 / 180.16

    def insulin_intake(self, t):
        """
        Insulin infusion. Default 1U at 30 minutes during 2 minutes.

        > t: time [min].
        < u: insulin intake [mU / min].
        """
        if self.insulin_intake_fun == None:    
            # Infusion of 1U during 5 minutes
            t0 = 30
            dt = 2
            u = 1000 / dt if (t >= t0 and t <= t0 + dt) else 0
            
            # Steady-state Basal Insulin
            u += self.u_basal
        else:
            u = self.insulin_intake_fun(t)
        return u
    
    def update_parameters(self, t):
        """
        Updates all oscillatory parameters of the subject.
        
        > t: time [min].
        < params: 8 parameters updated (only oscillatory ones).
        """
        new = [p_0 + np.sin(np.pi * t / (60 * U[0])) * p_0 * U[1] for p_0, U in zip(self.cyclic_params, self.circadian_time)]
        self.params[0:9] = new
        return new

    def full_model(self, t, state):
        """Scipy EDO 's solver function for the Hovorka 's Model"""
        # State Variables
        S1, S2, I, x1, x2, x3, G1, G2, Q1, Q2, C = state
        
        # Parameters
        _ = self.update_parameters(t)
        (EGP_0, F_01, k_12, S_id, S_ie, S_it, k_e, Bio, k_aint,
        k_a1, k_a2, k_a3, V_g, V_i, k_a, t_max, _, 
        R_th, R_cl, U_g_ceil) = self.params

        ### Auxiliar variables
        U_g = G2 / t_max
        t_max = G2 / U_g_ceil if U_g > U_g_ceil else t_max      # Time-to-maximum appearance rate of glucose
        G = Q1 / V_g                                            # Glucose
        EGP = max(EGP_0 * (1 - x3), 0)                          # Endogenous glucose production
        F_r = R_cl * (G - R_th) * V_g if G >= R_th else 0       # Renal glucose clearance
        #F_01_s = F_01 / 0.85                                    # 
        #F_01_c = F_01_s * G / (G + 1)                           # Non-insulin dependent glucose
        F_01_c = F_01 if G > 4.5 else F_01 * G / (4.5)           # Non-insulin dependent glucose

        # Food intake function
        D = self.food_intake(t) 
        
        # Insulin infusion function
        u = self.insulin_intake(t)

        # Debug
        #TEMP1 = -(F_01_c / (V_g * G) + x1) * Q1
        #TEMP2 = k_12 * Q2

        ### Differential equations of the model
        # Model order: Insulin Absorption -> Insulin Action, Food Absortion -> Glucose Dinamics -> Subcutaneous Glucose Dynamics
        derivs = np.array([ u - k_a * S1,                                                   # S1 - Masses os insulin in accessible compartment 
                            k_a * (S1 - S2),                                                # S2 - Masses os insulin in nonaccessible compartment
                            k_a * S2 / V_i - k_e * I,                                       # I - Plasma insulin concentration
                            -k_a1 * x1 + S_it * k_a1 * I,                                   # X1 - Distribution/transport effect
                            -k_a2 * x2 + S_id * k_a2 * I,                                   # X2 - Glucose disposal effect
                            -k_a3 * x3 + S_ie * k_a3 * I,                                   # X3 - Glucose production effect
                            -G1 / t_max + Bio * D,                                          # G1 - Glucose Masses in accessible compartment
                            (G1 - G2) / t_max,                                              # G2 - Glucose Masses in nonaccessible compartment
                            -(F_01_c / (V_g * G) + x1) * Q1 + k_12 * Q2 - F_r + EGP + U_g,  # Q1 - Glucose Masses in accessible compartment
                            x1 * Q1 - (k_12 + x2) * Q2,                                     # Q2 - Glucose Masses in nonaccessible compartment
                            k_aint * (G - C)], dtype='float64')                             # C - Glucose Subcutaneous Concentration in 
        # Return derivatives
        return derivs
    
    def steady_state(self, r = 7.0):
        """
        Solves the basal insulin for the Hovorka 's Model steady state.
        > r: basal glucose.
        < u: basal insulin.
        """
        # Parameters
        (EGP_0, F_01, k_12, S_id, S_ie, S_it, k_e, _, _,
        k_a1, k_a2, k_a3, V_g, V_i, _, _, _, 
        R_th, R_cl, _) = self.params
        
        ### Auxiliar variables
        F_r = R_cl * (r - R_th) * V_g if r >= R_th else 0       # Renal glucose clearance
        F_01_c = F_01 if r > 4.5 else F_01 * r / (4.5)           # Non-insulin dependent glucose
        k_b1 = S_it * k_a1
        k_b2 = S_id * k_a2
        k_b3 = S_ie * k_a3
        
        # Univariate function
        def q_u(u):
            return  (- F_01_c - F_r - r * V_g * k_b1 * u / (V_i * k_e * k_a1)
                     + k_12 * r * V_g * k_b1 * u / (k_a1 * V_i * k_e * (k_12 + k_b2 * u / (k_a2 * V_i * k_e)))
                     + EGP_0 * (1 - k_b3 * u / (V_i * k_e * k_a3)))
        
        # Univariate function derivative
        def dq_u(u):
            return  (- r * V_g * k_b1 / (V_i * k_e * k_a1)
                     + k_12 * r * V_g * k_b1 / (k_a1 * V_i * k_e * (k_12 + k_b2 * u / (k_a2 * V_i * k_e)))
                     + (-k_b2 * k_12 * r * V_g * k_b1 * u / (k_a2 * k_a1 * (V_i * k_e / k_12 + k_b2 / k_a2 * u) ** 2))
                     - EGP_0 * k_b3 / (V_i * k_e * k_a3))

        # Solve roots
        ss_root = root(q_u, 0, jac=dq_u, method='hybr')
        self.u_basal = ss_root.x

    #### Principais Métodos #####
    def simulate(self, y0, t0=0, tf=4*60, t_span=0.0005, cgm_reading_delay = 5):
        """
        Simulates the dynamics of the model for the current patient.

        > t0 (default 0): initial simulation time.
        > tf: final time to simulate (min).
        > t_span (default 0.0002): Time interval step.
        < solution: Array with the ODE solutions.
        < time: Array with time interval for each solution.
        """

        # Initial Conditions - S1, S2, I, x1, x2, x3, G1, G2, Q1, Q2, C
        basal_glucose = 7.0
        Q1 = self.params[12] * basal_glucose
        initial_conditions = np.array([0, 0, 6, 0, 0, 0, 0, 0, Q1, 0, basal_glucose])
        
        # Solving time interval
        t_calc = np.arange(t0, tf, t_span)

        # Solve ODE 's
        sol = solve_ivp(self.full_model,
                        (t0, tf),
                        initial_conditions,
                        method='RK45', t_eval=t_calc)
        print(sol.message)
        
        # Number of samples
        number_of_samples = int(tf * cgm_reading_delay)
        
        # Resamples the Continuous Glucose Monitor and its time
        glucose_split = np.array_split(sol.y[10, :], number_of_samples)
        time_split = np.array_split(sol.t, number_of_samples)
        self.cgm_glucose = np.array([i[0] for i in glucose_split])
        self.cgm_time =np.array([i[0] for i in time_split])
        
        # Adds Noise
        self.cgm_glucose = self.add_noise(self.cgm_glucose, type='ARMA')
        
        # Saves the result
        self.solution = sol.y
        self.time = sol.t

    def evaluate_metrics(self, tx = 10,ty = 4):
        """
        Returns a list of metrics evaluating the patient health based in his Blood Glucose levels.

        MBG: Mean Blood Glucose.
        SDBG: Standard Deviantion of Blood Glucose.
        TIZ: Time in zone.
        TA_tx: Time above threshold tx.
        TB_ty: Time bellow threshold ty.
        MAGE: Mean Amplitude of Glycemic Excursion TODO.
        ADRR: Average Daily Risk Range.
        ADRR_risk: Risk according to ADRR.
        LBGI: Low Blood Glucose Index.
        HBGI: High Blood Glucose Index.
        HbA1c: Glycated hemoglobin.

        < dictionary with the metric described above.
        """
        # Blood Glucose Results
        bg = self.solution[8, :] / self.params[12]

        # Time span between values
        time_span = self.time[1] - self.time[0]
        
        # TIZ, TA_tx, TA_ty
        TA_tx = len(bg[bg > tx]) * time_span
        TB_ty = len(bg[bg < ty]) * time_span
        TIZ = len(bg) * time_span - TA_tx - TB_ty

        # MAGE
        MAGE = None

        # ADRR, LBGI, HBGI
        f_bg = 1.509 * (np.log(bg * 18.0182)**1.084 - 5.381)
        rl = np.array([10 * (f_bg_i ** 2) if f_bg_i < 0 else 0 for f_bg_i in f_bg])
        rh = np.array([10 * (f_bg_i ** 2) if f_bg_i > 0 else 0 for f_bg_i in f_bg])
        LR_l = np.max(rl)
        HR_l = np.max(rh)
        ADRR = (LR_l + HR_l) / 2.0      # TODO convert to each day
        LBGI = np.mean(rl)
        HBGI = np.mean(rh)
        if ADRR < 20:
            ADRR_risk = 'Low'
        elif ADRR > 40:
            ADRR_risk = 'High'
        else:
            ADRR_risk = 'Medium'

        # MBG, SDBG
        MBG = np.mean(bg)
        SDBG = np.std(bg)

        # HbA1c
        HbA1c_1 = (MBG + 4.29) / 1.98
        HbA1c_2 = (MBG + .1) / 1.54
        HbA1c_3 = (MBG - .47) / 1.23
        HbA1c_4 = (MBG + 3.61) / 1.87

        HbA1c = np.mean([HbA1c_1, HbA1c_2, HbA1c_3, HbA1c_4])

        metrics = { 'MBG': MBG, 'SDBG': SDBG, 'TIZ': TIZ,
                    'TA_tx': TA_tx, 'TB_ty': TB_ty, 'MAGE': MAGE,
                    'ADRR': ADRR, 'ADRR_Risk': ADRR_risk,'LBGI': LBGI, 'HBGI': HBGI, 
                    'HbA1c': HbA1c, 'Risk_Space':f_bg}
        
        return metrics

    @staticmethod
    def add_noise(array, type='GAUSS'):
        """
        Adds noise in given a array.
        > array: array in which the noise will be added.
        > type: type of noise ('GAUSS' ou 'ARMA')
        
        <array: noisy array.
        """
        if type == 'ARMA':
            
            # ARMA noise parameters
            size = array.shape[0]
            e_n = np.empty(size)
            e_n[0] = np.random.normal()
            lamb = 15.96
            eps = -5.471
            delta = 1.6898
            gamma = -0.5444
            for i in range(size - 1):
                e_n[i + 1] = 0.7 * (e_n[i] + np.random.normal())
            
            error_n = eps + lamb * np.sinh((e_n - gamma) / delta)
            assert(size == error_n.shape[0])

            return array + error_n * 0.0555 # converts to mmol/L
        else:
            return array + np.random.normal()

    def plot_results(self, fig=1, insulin=False, sc_glucose=True, cgm_glucose=True, bolus_insulin=True, meal=True,
                        aux=[False, False, False, False, False], all=False, parameteres=False, total_samples=60):
        """
        Plota os gráficos da Simulação, de acordo com a entrada.

        -insulin: Se o gráfico da insulina será plotado.
        -glucose: Se o gráfico da glicose será plotado.
        -sc_glucose: Se o gráfico da glicose subcutânea será plotado.
        -aux: lista com os gráficos auxiliares que serão postados. 
        A ordem das variáveis são: S1, S2, I, x1, x2, x3, G1, G2, Q1, Q2, C.
        -all: se todos os gráficos serão plotados.
        """
        # Axes count
        axes_c = int(insulin) + 1

        # Convertstime to TimeStamps
        datetime_list = [datetime.datetime(2018, 8, 20, int(t / 60), int(t % 60)) for t in self.time.squeeze()]
        plot_time = md.date2num(datetime_list)
        
        # CGM glucose data
        glucose_cgm = self.cgm_glucose
        time_cgm = self.cgm_time
        # CGM Glucose date
        date_cgm = [datetime.datetime(2018, 8, 20, int(t / 60), int(t % 60)) for t in time_cgm]

        # Formatting Timestamps
        formatter = md.DateFormatter('%H:%M')
        locator = md.HourLocator(interval = 1)

        # Inicialize figures
        fig_main = plt.figure(fig, figsize=(8, 6)) # Para o plot da insulina, glucose/sc_glucose
        if all or np.any(aux) or parameteres:
            fig_aux = plt.figure(fig + 1)              # Para as demais funções
        
        # Eixos da insulina e glicose
        G_ax = fig_main.add_subplot(axes_c, 1, 1)
        if insulin:
            I_ax = fig_main.add_subplot(axes_c, 1, 2, sharex=G_ax)
        
        # Eixo dos parâmetros
        if parameteres:
            Param_ax = fig_aux.add_subplot(111)
            Param_ax.xaxis.set_major_formatter(formatter)
            Param_ax.xaxis.set_major_locator(locator)
        # Principais plots
        # Glicose
        G_ax.set_xlabel('Tempo [H:M]')
        G_ax.set_ylabel('Concentração de Glicose [mmol/L]')
        #G_ax.set_ylim((4, 10))
        G_ax.xaxis.set_major_formatter(formatter)
        G_ax.xaxis.set_major_locator(locator)
        
        # Glicose do plasma
        _, = G_ax.plot(plot_time, self.solution[8,:] / self.params[12], color='blue', label="Glicose no Plasma")
        
        # Glicose subcutânea
        if sc_glucose:
            _, = G_ax.plot(plot_time, self.solution[10,:], color='red', label="Glicose Subcutânea")    
        
        # Glicose do cgm (subcutânea com ruído)
        if cgm_glucose:
            y = glucose_cgm
            G_ax.scatter(date_cgm, y, marker='.', color='green', label="Leitura do Sensor")

        # Insulina
        if insulin:
            I_ax.plot(plot_time, self.solution[2, :], color='blue')
            I_ax.set_xlabel('Horas [min]')
            I_ax.set_ylabel('Concentração de Insulina [mU/L]')

        # Parâmetros15:24 06/11/2018
        if parameteres:
            # Exemplo TODO
            parameters_y = [self.update_parameters(t)[0] for t in self.time]
            Param_ax.plot(plot_time, parameters_y, color='blue')
            #I_ax.set_xlabel('Horas [min]')
            Param_ax.set_ylabel('Valor do Parâmetro')
        
        plt.legend(loc='best')
        plt.draw()

# Main call
S1 = Subject()
S1.simulate(0, tf = 60 * 12)
S1.plot_results(1)
print(S1.u_basal)
print(S1.evaluate_metrics())
# S2 = Subject()
# S2.simulate(0, tf=60 * 12)
# S2.plot_results(2)

# S3 = Subject()
# S3.simulate(0, tf=60 * 12)
# S3.plot_results(3)

plt.show()
