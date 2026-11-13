import numpy as np
from pymoo.core.problem import ElementwiseProblem
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize
from pymoo.visualization.scatter import Scatter
import matplotlib.pyplot as plt

RHO = 7850.0                                    # steel, kg/m^3 (swap for your material)
BETA1, BETA2 = np.radians(20), np.radians(160)  # Gimbal angle range for reference gimbal speed
TAU_HOLD, T_HOLD = 0.7, 0.4                     # Hold a torque of >= 0.7 Nm for 0.4 s
BETA_DOT_REF = (BETA2-BETA1)/T_HOLD             # Gimbal speed to hold torque for 0.4 s
SIGMA_ALLOW = 200e6                             # Allowable stress of steel, Pa
KE_MAX_WHEEL = 100                              # Maximum kinetic energy of the wheel, J
KE_AREA_MAX = np.inf                            # Maximum kinetic energy per unit area of the wheel, J/m^2
TAU_PEAK_MAX = 3.0                              # Maximum peak torque of the scissor pair, N·m
T_RIM_MIN, B_MAX, R_MAX = 0.002, 0.030, 0.04    # Minimum rim thickness, maximum web width, maximum outer radius, m
BETA_DDOT_MAX = 5000                            # Maximum gimbal angular acceleration, rad/s^2
THETA_DOT_MAX = 20000*np.pi/30                  # Maximum flywheel angular velocity, rad/s
OMEGA_LIMB = 5.0                                # Limb angular velocity, rad/s
T_GIMBAL = 0.35                                 # Gimbal torque, N·m
BETA_DOT_MAX = 350*np.pi/30                     # Max Gimbal motor speed, rad/s
J_WHEEL_MOTOR = 0.09*0.005**2                   # Flywheel motor inertia in gimbal axis, kg·m^2
J_GIMBAL_MOTOR = 61*10**-7                      # Gimbal motor inertia, kg·m^2
NU = 0.3

def wheel_props(x):
    R, r, B, b, theta_dot, Beta_dot, n, Beta_ddot = x
    m_rim = np.pi*RHO*B*(R**2 - r**2)
    m_web = np.pi*RHO*b*r**2
    J = 0.5*np.pi*RHO*B*(R**4 - r**4) + 0.5*np.pi*RHO*b*r**4
    sigma = RHO*theta_dot**2/4*((1 - NU)*r**2 + (3 + NU)*R**2)   # Pa, peak tangential stress
    return {
        "mass_wheel": m_rim + m_web,          # kg, one wheel
        "J": J,                               # kg·m², about the spin axis, one wheel
        "H": J*theta_dot,                     # N·m·s, angular momentum per wheel
        "tip_speed": R*theta_dot,             # m/s
        "sigma": sigma,                       # Pa
        "sigma_ring": RHO*(R*theta_dot)**2,   # Pa, thin-ring upper bound (ρv²)
    }

class CMG(ElementwiseProblem):
    def __init__(self):
        
        super().__init__(n_var=8, n_obj=2, n_ieq_constr=10,
                         xl=[0.010, 0.0015, 0.005, 0.005, 200, 0.5, 1, 0.1],
                         xu=[R_MAX, R_MAX, B_MAX, B_MAX, THETA_DOT_MAX, 150, 5, BETA_DDOT_MAX])
        #x= [Outer radius, rim thickness, rim width, web width, flywheel velocity, gimbal velocity, gear ratio, gimbal acceleration]
        
    def _evaluate(self, x, out, *args, **kwargs):
        R, r, B, b, theta_dot, Beta_dot, n, Beta_ddot = x

        p = wheel_props(x)
        mass, J, sigma = p["mass_wheel"], p["J"], p["sigma"]
        #mass = np.pi*RHO*B*(R**2-r**2) + np.pi*RHO*b*r**2               # mass of the scissor pair, kg
        #J = 0.5*np.pi*RHO*B*(R**4-r**4) + 0.5*np.pi*RHO*b*r**4          # inertia of the scissor pair about the gimbal axis, kg·m^2
        tau_peak = 2 * J* theta_dot * Beta_dot                          # peak torque of the scissor pair at 90 degrees, N·m
        #sigma = RHO*theta_dot**2/4*((1 - 0.3)*r**2 + (3 + 0.3)*R**2)    # peak stress in the flywheel rim, Pa
        KE_wheel = 1/4*J*(2*theta_dot**2 + Beta_dot**2)                 # kinetic energy of the scissor pair, J
        KE_area = KE_wheel/(2*np.pi*R*(R+B))                            # kinetic energy per unit area of the scissor pair, J/m^2
        J_min = TAU_HOLD/(2*BETA_DOT_REF*theta_dot*np.sin(BETA1))       # Minimum inertia of flywheel for TAU_HOLD to last T_HOLD seconds at reference gimbal speed
        T_mg = Beta_ddot*(1.5*(0.5*J+J_WHEEL_MOTOR)/n + J_GIMBAL_MOTOR*n) + J*theta_dot*OMEGA_LIMB/n    # Gimbal motor torque 
        Beta_dot_mg = n*Beta_dot                                        # Gimbal motor speed 
        
        # g <= 0 constraints
        g_1 = 1 - J/J_min
        g_2 = KE_wheel/KE_MAX_WHEEL - 1
        g_3 = KE_area/KE_AREA_MAX - 1
        g_4 = sigma/SIGMA_ALLOW -1
        g_5 = tau_peak/TAU_PEAK_MAX - 1
        g_6 = (b - B)/B_MAX
        g_7 = (T_RIM_MIN - (R - r))/R_MAX
        g_8 = Beta_dot**2/(np.pi*Beta_ddot) - 1
        g_9 = T_mg/T_GIMBAL - 1
        g_10 = Beta_dot_mg/BETA_DOT_MAX - 1

        g = np.array([g_1, g_2, g_3, g_4, g_5, g_6, g_7, g_8, g_9, g_10])
        
        out["F"] = [2 * mass, -tau_peak]                          # minimize mass, minimize KE
        out["G"] = g

res = minimize(CMG(), NSGA2(pop_size=200), ("n_gen", 300), seed=1, verbose=False)
if res.F is None:
    raise SystemExit("No feasible design found. Loosen KE_MAX_WHEEL or another limit.")

problem = CMG()
F = res.F.copy()
F[:, 1] = -F[:, 1]                      # un-negate torque for plotting
names = ["R (m)", "r (m)", "B (m)", "b (m)", "theta_dot (rad/s)", "beta_dot (rad/s)", "n (unitless)", "beta_ddot (rad/s^2)"]

fig, ax = plt.subplots()
ax.scatter(F[:, 0], F[:, 1], picker=True, pickradius=5)
ax.set_xlabel("Total flywheel mass (kg)")
ax.set_ylabel("Peak torque (N·m)")
ax.grid(True)

def on_pick(event):
    for i in event.ind:
        out = problem.evaluate(res.X[i], return_as_dictionary=True)
        p = wheel_props(res.X[i])
        print(f"\nPoint {i}: mass={F[i,0]:.4f} kg, torque={F[i,1]:.3f} N·m")
        for n, v in zip(names, res.X[i]):
            print(f"  {n:20s} {v:.5g}")
        print(f"  {'J (per wheel)':20s} {p['J']:.4e} kg·m²")
        print(f"  {'H (per wheel)':20s} {p['H']:.4f} N·m·s")
        print(f"  {'tip speed':20s} {p['tip_speed']:.1f} m/s")
        print(f"  {'stress (hollow disk)':20s} {p['sigma']/1e6:.1f} MPa  (limit {SIGMA_ALLOW/1e6:.0f})")
        print(f"  {'stress (thin ring)':20s} {p['sigma_ring']/1e6:.1f} MPa")
        print(f"  {'safety margin':20s} {SIGMA_ALLOW/p['sigma']:.2f}x")
        print("  constraints G (<=0 ok, ~0 active):", np.round(out["G"], 4))

fig.canvas.mpl_connect("pick_event", on_pick)
plt.show()
