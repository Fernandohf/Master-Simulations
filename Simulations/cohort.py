"""
File with the main Classes.
"""
import matplotlib.pyplot as plt
import matplotlib.dates as md
import matplotlib.collections as collections
from scipy.integrate import solve_ivp
from scipy.optimize import root
import numpy as np
import datetime

# Hovorka 's Model ########################
# TODO LIST
# - TRANSLATE
# - Write Tests
# - CONTROL
# - Unknown Error makes simulation fail - likely parameters


class Subject():
    """
    Subject Simulation Class.
    """
    # Paramete Distributions
    # Paper Software for in Silico Testing.pdf
    PARAM_DIST = {"EGP_0": {"mean": 0.0169, "std": 0.0039, "unit": "mmol/min", "distribution": "Gaussian", "variability": [[0, 3], 0.05]},
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
                  "V_g": {"mean": -1.89711998489, "std": 0.23, "unit": "L/kg", "distribution": "Lognormal", "variability": None},
                  "V_i": {"mean": 0.12, "std": 0.012, "unit": "L/kg", "distribution": "Gaussian", "variability": None},
                  "k_a": {"mean": 0.018, "std": 0.0045, "unit": "mmol/L", "distribution": "Gaussian", "variability": [[0, 3], 0.05]},
                  "t_max": {"mean": -3.689, "std": 0.025, "unit": "min", "distribution": "Lognormal", "variability": None},
                  "w": {"mean": 74.9, "std": 14.4, "unit": "kg", "distribution": "Gaussian", "variability": None},
                  "R_th": {"mean": 9, "std": 1.5, "unit": "mmol/L", "distribution": "Gaussian", "variability": None},
                  "R_cl": {"mean": 0.01, "std": 0.025, "unit": "1/min", "distribution": "Gaussian", "variability": None},
                  "U_ceil": {"mean": 0.02, "std": 0.035, "unit": "TODO ??", "distribution": "Uniform", "variability": None}
                  }

    def __init__(self, params=None, food_intake_f=None, control_type=None):
        """
        Initializes the Class Subject.

        > params(default None): list with all the parameters in Hovorka 's Model in order given bellow.
        In case of None, automatic parameters will be given following the built-in probability distribution.
        [EGP_0, F_01, k_12, S_id, S_ie, S_it, k_e, Bio, k_a1, k_a2, k_a3, k_a_int, V_g, V_i, k_a, t_max, w, R_th, R_cl, U_ceil]

        > food_intake_f(default None): Function in format 'fun(t)' that returns the rate of carbohydrate ingestetion in given time t.
        The values should be (g / min).

        > insulin_intake_f(default None): Function in format 'fun(t)' that return the rate of  insulin infusion in given time t.
        The values should be (U / min).

        """
        # Constants
        self.CGM_AVERAGE_DELAY = 15  # minutes
        self.SAMPLING_TIME = 15      # minutes

        # Last Control Variables
        self.control_u = np.array([0])          # output after saturation
        self.control_I = np.array([0])          # previous I
        self.control_v = np.array([0])          # output before saturation

        # Variable with resulting simulation values from 'simulate'
        self.solution = np.array([])
        self.time = np.array([])
        self.max_glucose = 0
        self.min_glucose = 0

        # Steady-state basal insulin value
        self.u_basal = 0

        # Auxiliary input functions: food intake and control infusion
        self.control_type = control_type
        self.food_intake_fun = food_intake_f

        # Continuous Variables
        self.insulin = np.array([])
        self.meal = np.array([])

        # Sensor reading
        self.cgm_glucose = np.array([])
        self.cgm_time = np.array([])

        # Parameters Initialization
        if params is None:
            # Initalize parameters from defined distributions
            # [   *0,   *1,   *2,   *3,   *4,   *5,  *6,  *7,     *8,    9,   10,  11,   12,  13,  14,    15,16,   17,   18,     19]
            # [EGP_0, F_01, k_12, S_id, S_ie, S_it, k_e, Bio, k_aint, k_a1, k_a2, k_a3, V_g, V_i, k_a, t_max, w, R_th, R_cl, U_ceil]
            self.params = self.get_params_from_dist(self.__class__.PARAM_DIST)
        else:
            self.params = params

        # Steady-state basal insulin
        self.steady_state()

    def get_params_from_dist(self, param_dist_dict):
        """
        Collects the values from the parameters distributions.

        The dictionary should be a dictionary of dictionaries and should have all these fields:

        "EGp_dict, F_01, k_12, S_id, S_ie, S_it, k_e, Bio, k_aint, k_a1, k_a2, k_a3, V_g, V_i, k_a, t_max, w, R_th, R_cl, U_ceil"

        And each field should result in a dictionary with the fields bellow:
        "mean": 0.0169, "std": 0.0039, "unit": "mmol/min", "distribution": "Gaussian", "variability": "Oscillatory"
        Description:
        - mean": Average values of the distribution. OBS.: If the distribution is uniform, this is the min value.
        - std": Standard deviation for the distribution. OBS.: If the distribution is uniform, this is the max value.
        - unit": Units of the parameter.
        - distribution": The type of distribution, ex.: "Gaussian", "Lognormal", etc.
        - variability: List in the format [ min oscillatory frequency, max oscillatory frequency], max amplitude ] if the parameter is
            oscillatory, None otherwise.

        The dictionary returned has all these fields:

        "EGp_dict, F_01, k_12, S_id, S_ie, S_it, k_e, Bio, k_aint, k_a1, k_a2, k_a3, V_g, V_i, k_a, t_max, w, R_th, R_cl, U_ceil"

        And each inner dictionary should look like this:
        "value": 0.0169, "var_freq": 2, "var_amp": 0.05, "unit": "mmol/min"
        Description:
        - value: Average values of the distribution.
        - var_freq: Sinusoidal wavelength (frequency) of variability.
        - var_amp: Sinusoidal amplitude of variability.

        > param_dist_dict: Dictionary with all the parameter distributions.

        < param_values_dict: Dictionary with all the values from the distributions.

        """
        param_values_dict = {}
        # For each item in the dict
        for params, dist_dict in param_dist_dict.items():

            # Initialize inner dict
            param_values_dict[params] = {}

            # For each type of distribution
            dist_type = dist_dict["distribution"]
            if dist_type == "Gaussian" or dist_type == "Normal":
                param_values_dict[params]["value"] = np.random.normal(dist_dict["mean"], dist_dict["std"])
            elif dist_type == "Lognormal":
                param_values_dict[params]["value"] = np.random.lognormal(dist_dict["mean"], dist_dict["std"])
            elif dist_type == "Uniform":
                param_values_dict[params]["value"] = np.random.uniform(dist_dict["mean"], dist_dict["std"])

            # For each Stationary and Oscillatory variability
            var_type = dist_dict["variability"]
            if var_type:
                # First value of variability min
                param_values_dict[params]["var_freq"] = np.random.uniform(var_type[0][0], var_type[0][1])
                param_values_dict[params]["var_amp"] = var_type[1]
            else:
                param_values_dict[params]["var_freq"] = 1
                param_values_dict[params]["var_amp"] = 0.0

        # Caveats
        w = param_values_dict["w"]["value"]
        param_values_dict["t_max"]["value"] = 1. / param_values_dict["t_max"]["value"]
        param_values_dict["V_g"]["value"] = param_values_dict["V_g"]["value"] * w
        param_values_dict["V_i"]["value"] = param_values_dict["V_i"]["value"] * w
        param_values_dict["U_ceil"]["value"] = param_values_dict["U_ceil"]["value"] * w

        # Return final values
        return param_values_dict

    def update_and_unpack_params(self, t):
        """
        Unpacks parameters and updates oscillatory ones according to value of t.
        EGP_0, F_01, k_12,  S_id, S_ie, S_it,  k_e, Bio,  k_aint, k_a1,
        k_a2,  k_a3, V_g,   V_i,  k_a,  t_max, w,   R_th, R_cl,   U_ceil

        < values_list: List of the unpacked values.
        """
        # For each parameter
        parameters = ["EGP_0", "F_01", "k_12", "S_id", "S_ie", "S_it", "k_e", "Bio", "k_a_int", "k_a1",
                      "k_a2", "k_a3", "V_g", "V_i", "k_a", "t_max", "w", "R_th", "R_cl", "U_ceil"]
        values_list = []
        for p in parameters:
            p_dict = self.params[p]
            values_list.append(p_dict["value"] + np.sin(np.pi * t / (60 * p_dict["var_freq"])) * p_dict["value"] * p_dict["var_amp"])

        # Returns values in given order
        return tuple(values_list)

    def food_intake(self, t):
        """
        Food ingestion function. Default 50g at 30 minutes during 15 minutes.

        > t: time [min].
        < D: Equivalent glicose ingested [mmol /min].
        """

        if self.food_intake_fun is None:
            if t > 30 and t < 45:
                # 50g during 15 minutes
                t0 = 30
                dt = 15
                d = 50 / dt
            elif t > 260 and t < 280:
                # 100g during 15 minutes
                t0 = 260
                dt = 20
                d = 100 / dt
            elif t > 520 and t < 540:
                # 70g during 15 minutes
                t0 = 520
                dt = 20
                d = 70 / dt
            else:
                d = 0
        else:
            d = self.food_intake_fun(t)

        # Convert from [g/min] to [mmol/min]
        return d * 1000 / 180.16

    def insulin_intake(self, t, C):
        """
        Insulin infusion. Default 1U at 30 minutes during 2 minutes.

        > t: time [min].
        < u: insulin intake [mU / min].
        """
        if self.control_type is None:
            # Infusion of 1U during 5 minutes
            t0 = 30
            dt = 5
            u = 1000 / dt if (t >= t0 and t <= t0 + dt) else 0

            # Steady-state Basal Insulin
            u = self.u_basal
        elif self.control_type == 'PID':
            u = self.controller_PID(C)
        return u

    def full_model(self, t, state):
        """Scipy EDO 's solver function for the Hovorka 's Model"""
        # State Variables
        S1, S2, I, x1, x2, x3, G1, G2, Q1, Q2, C = state

        # Updates Oscillatory and get parameters values
        (EGP_0, F_01, k_12, S_id, S_ie, S_it,
         k_e, Bio, k_aint, k_a1, k_a2, k_a3,
         V_g, V_i, k_a, t_max, _, R_th, R_cl,
         U_ceil) = self.update_and_unpack_params(t)

        # Auxiliar variables
        U_g = G2 / t_max
        t_max = G2 / U_ceil if U_g > U_ceil else t_max      # Time-to-maximum appearance rate of glucose
        G = Q1 / V_g                                        # Glucose
        EGP = max(EGP_0 * (1 - x3), 0)                      # Endogenous glucose production
        F_r = R_cl * (G - R_th) * V_g if G >= R_th else 0   # Renal glucose clearance
        F_01_s = F_01 / 0.85
        F_01_c = F_01_s * G / (G + 1)                       # Non-insulin dependent glucose
        # F_01_c = F_01 if G > 4.5 else F_01 * G / (4.5)     # Non-insulin dependent glucose

        # Food intake function
        D = self.food_intake(t)

        # Insulin infusion function
        u = self.insulin_intake(t, C)

        # Differential equations of the model
        # Model order: Insulin Absorption -> Insulin Action, Food Absortion -> Glucose Dinamics -> Subcutaneous Glucose Dynamics
        derivs = np.array([u - k_a * S1,                                                   # S1 - Masses os insulin in accessible compartment
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

    def steady_state(self, r=5.0):
        """
        Solves the basal insulin for the Hovorka 's Model steady state.
        > r: basal glucose.
        < u: basal insulin.
        """
        # Parameters
        (EGP_0, F_01, k_12, S_id, S_ie, S_it, k_e, _, _, k_a1, k_a2, k_a3,
         V_g, V_i, _, _, _, R_th, R_cl, _) = self.update_and_unpack_params(0)

        # Auxiliar variables
        F_r = R_cl * (r - R_th) * V_g if r >= R_th else 0       # Renal glucose clearance
        F_01_c = F_01 if r > 4.5 else F_01 * r / (4.5)          # Non-insulin dependent glucose
        k_b1 = S_it * k_a1
        k_b2 = S_id * k_a2
        k_b3 = S_ie * k_a3

        # Univariate function
        def q_u(u):
            return (-F_01_c - F_r - r * V_g * k_b1 * u / (V_i * k_e * k_a1) +
                    k_12 * r * V_g * k_b1 * u / (k_a1 * V_i * k_e * (k_12 + k_b2 * u / (k_a2 * V_i * k_e))) +
                    EGP_0 * (1 - k_b3 * u / (V_i * k_e * k_a3)))

        # Univariate function derivative
        def dq_u(u):
            return (-r * V_g * k_b1 / (V_i * k_e * k_a1) +
                    k_12 * r * V_g * k_b1 / (k_a1 * V_i * k_e * (k_12 + k_b2 * u / (k_a2 * V_i * k_e))) +
                    (-k_b2 * k_12 * r * V_g * k_b1 * u / (k_a2 * k_a1 * (V_i * k_e / k_12 + k_b2 / k_a2 * u) ** 2)) -
                    EGP_0 * k_b3 / (V_i * k_e * k_a3))

        # Solve roots
        ss_root = root(q_u, 0, jac=dq_u, method='hybr')
        self.u_basal = ss_root.x

    def partial_simulation(self, t0, tf, initial_conditions, t_step=0.002):

        # Solving time interval
        t_calc = np.arange(t0, tf, t_step)

        # Solve ODE 's
        sol = solve_ivp(self.full_model,
                        (t0, tf),
                        initial_conditions,
                        method='RK45', t_eval=t_calc)
        # Check if the solution was successfully found
        if sol.status == -1:
            print(sol.message)
        return sol

    def simulate(self, y0, t0=0, tf=4 * 60, t_span=5):
        """
        Simulates the dynamics of the model for the current patient.

        > t0 (default 0): initial simulation time.
        > tf: final time to simulate (min).
        > t_span (default 5): Time interval between steps.
        < solution: Array with the ODE solutions.
        < time: Array with time interval for each solution.
        """

        # Initial Conditions - S1, S2, I, x1, x2, x3, G1, G2, Q1, Q2, C
        initial_glucose = 5.
        Q1 = self.params["V_g"]["value"] * initial_glucose
        initial_conditions = np.array([0, 0, self.u_basal, 0, 0, 0, 0, 0, Q1, 0, initial_glucose])

        # Placegolder
        solutions = initial_conditions.copy().reshape(-1, 1)
        time_sol = np.array(t0)
        cgm_reading = np.array([initial_glucose])
        cgm_time = np.array([t0])

        # For each t_span minutes
        for t in range(t0, tf - t_span, t_span):
            solution = self.partial_simulation(t, t + t_span, initial_conditions)

            # Last values from the solver
            initial_conditions = solution.y[:, -1]

            # Append results
            solutions = np.hstack((solutions, solution.y))
            time_sol = np.hstack((time_sol, solution.t))
            cgm_reading = np.hstack((cgm_reading, solution.y[-1, -1]))
            cgm_time = np.hstack((cgm_time, solution.t[-1]))

        # Saves results
        self.solution = solutions
        self.cgm_glucose = cgm_reading
        self.cgm_time = cgm_time
        self.time = time_sol

        # Adds Noise to CGM readings
        self.cgm_glucose = self.add_noise(self.cgm_glucose, type='ARMA')

    @staticmethod
    def solve_metric_LHBGI(bg):
        """
        TODO A method to solve each metrics.

        """
        f_bg = 1.509 * (np.log(bg * 18.0182)**1.084 - 5.381)
        rl = np.array([10 * (f_bg_i ** 2) if f_bg_i < 0 else 0 for f_bg_i in f_bg])
        rh = np.array([10 * (f_bg_i ** 2) if f_bg_i > 0 else 0 for f_bg_i in f_bg])
        LR_l = np.max(rl)
        HR_l = np.max(rh)

        LBGI = np.mean(rl)
        HBGI = np.mean(rh)
        return LBGI, HBGI

    def evaluate_metrics(self, tx=10, ty=4):
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
        bg = self.solution[8, :] / self.params["V_g"]["value"]

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
        ADRR = (LR_l + HR_l) / 2.0  # TODO convert to each day
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

        metrics = {'MBG': MBG, 'SDBG': SDBG, 'TIZ': TIZ,
                   'TA_tx': TA_tx, 'TB_ty': TB_ty, 'MAGE': MAGE,
                   'ADRR': ADRR, 'ADRR_Risk': ADRR_risk, 'LBGI': LBGI, 'HBGI': HBGI,
                   'HbA1c': HbA1c, 'Risk_Space': f_bg}

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

            return array + error_n * 0.0555  # converts to mmol/L
        else:
            return array + np.random.normal()

    @staticmethod
    def saturation(input):
        """
        Limits the input for the given pump range: 0 < input < 300 mU.

        > input: value.
        < output: restricted value.
        """
        return max(0, min(input, 300))

    def controller_PID(self, y, d=0):
        """
        Literature PI Controller. Software paper
        > y: Process variable.
        > d: Optional feedforward metric

        < u: Control signal.
        """
        # Target
        r = 5
        # Proportional term
        K_p = -20
        K_i = -55
        K_d = -2

        P = K_p * (r - y)

        # Integrative term
        I = (self.control_I[-1] + K_i * (r - y) / self.SAMPLING_TIME +
             1 / self.SAMPLING_TIME * (self.control_u[-1] - self.control_v[-1]))

        # Feed Forward
        CR = .8
        FF = K_d * CR * d

        # Output
        v = P + I + FF
        u = self.saturation(v)

        # Update values
        np.append(self.control_I, I)
        np.append(self.control_v, v)
        np.append(self.control_u, u)

        # Safety Logic
        if y < 4:
            u = 0

        return u

    def controller_SM(self):
        pass

    def controller_MPC(self):
        pass
    # TODO Add controller class? maybe

    def plot_results(self, fig=1, insulin=False, sc_glucose=True, cgm_glucose=True, bolus_insulin=True, meal=True,
                     aux=[False, False, False, False, False], all=False, metrics=True, total_samples=60):
        """
        Plota os gráficos da Simulação, de acordo com a entrada.

        -insulin: Se o gráfico da insulina será plotado.
        -glucose: Se o gráfico da glicose será plotado.
        -sc_glucose: Se o gráfico da glicose subcutânea será plotado.
        -aux: lista com os gráficos auxiliares que serão postados.
        A ordem das variáveis são: S1, S2, I, x1, x2, x3, G1, G2, Q1, Q2, C.
        -all: se todos os gráficos serão plotados.
        """
        # Constants
        INTERVAL = 4  # hours

        # Axes count
        axes_c = int(insulin) + 1 + int(metrics)

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
        locator = md.HourLocator(interval=INTERVAL)

        # Inicialize figures
        fig_main = plt.figure(fig, figsize=(8, 6))     # Para o plot da insulina, glucose/sc_glucose
        if all or np.any(aux) or metrics:
            fig_aux = plt.figure(fig + 1)              # Para as demais funções

        # Eixos da insulina e glicose
        G_ax = fig_main.add_subplot(axes_c, 1, 1)
        if insulin:
            I_ax = fig_main.add_subplot(axes_c, 1, 2, sharex=G_ax)

        # Eixo dos parâmetros
        if metrics:
            metrics_ax = fig_main.add_subplot(axes_c, 1, 2, sharex=G_ax)

        # Principais plots

        # Glicose
        G_ax.set_xlabel('Tempo [H:M]')
        G_ax.set_ylabel('Concentração de Glicose [mmol/L]')
        #G_ax.set_ylim((4, 10))
        G_ax.xaxis.set_major_formatter(formatter)
        G_ax.xaxis.set_major_locator(locator)

        # Glicose do plasma
        bg = self.solution[8, :] / self.params["V_g"]["value"]
        _, = G_ax.plot(plot_time, bg, color='blue', label="Glicose no Plasma")

        # Glicose subcutânea
        if sc_glucose:
            _, = G_ax.plot(plot_time, self.solution[10, :], color='red', label="Glicose Subcutânea")

        # Glicose do cgm (subcutânea com ruído)
        if cgm_glucose:
            y = glucose_cgm
            G_ax.scatter(date_cgm, y, marker='.', color='green', label="Leitura do Sensor")

        # Insulina
        if insulin:
            # Insulin Concentration
            I_ax.plot(plot_time, self.solution[2, :], color='blue')
            I_ax.set_xlabel('Tempo [H:M]')
            I_ax.set_ylabel('Concentração de Insulina [mU/L]')

            # Target insulin
            I_ax.plot(plot_time, np.repeat(self.u_basal, len(plot_time)), color='lime')

        # HLBGI Metrics
        div = 80
        n_splits = (len(bg) // div)
        remaining = len(bg) - (n_splits * div)
        bg_sliced = np.split(bg[:n_splits * div], div)

        L_data = np.array([])
        H_data = np.array([])
        for bg_step in bg_sliced:
            L_metric, H_metric = self.solve_metric_LHBGI(bg_step)
            L_data = np.concatenate((L_data, np.repeat(L_metric, len(bg_step))))
            H_data = np.concatenate((H_data, np.repeat(H_metric, len(bg_step))))
        L_last = L_data[-1]
        H_last = H_data[-1]
        L_data = np.concatenate((L_data, np.repeat(L_last, remaining)))
        H_data = np.concatenate((H_data, np.repeat(H_last, remaining)))

        metrics_ax.plot(plot_time, -L_data, color='red')
        metrics_ax.plot(plot_time, H_data, color='purple')
        metrics_ax.set_xlabel('Tempo [H:M]')
        metrics_ax.set_ylabel('Blood Glucose Index')
        # Colletions
        # collection = collections.BrokenBarHCollection.span_where(
        #     plot_time, ymin=0, ymax=max(H_data), where=H_data > 0, facecolor='blue', alpha=0.3)
        # metrics_ax.add_collection(collection)
        # collection = collections.BrokenBarHCollection.span_where(
        #     plot_time, ymin=min(L_data), ymax=0, where=L_data < 0, facecolor='red', alpha=0.3)
        # metrics_ax.add_collection(collection)

        # Legends
        G_ax.legend()
        # I_ax.legend()

        plt.legend(loc='best')
        plt.draw()

# Main call
np.random.seed(12)
S1 = Subject(control_type='PID')
S1.simulate(0, tf=60 * 24)
S1.plot_results(1, insulin=False, metrics=True)
print(S1.u_basal)
print(S1.evaluate_metrics())
# S2 = Subject()
# S2.simulate(0, tf=60 * 12)
# S2.plot_results(2)

# S3 = Subject()
# S3.simulate(0, tf=60 * 12)
# S3.plot_results(3)

plt.show()
