import numpy as np
import matplotlib.pyplot as plt
import math

class gyro:

    def __init__(self, gimbalspeed, rotospeed, rho, R, r, r_s, B, b):
        self.Beta_dot = gimbalspeed
        self.Theta_dot = rotospeed
        self.rho = rho
        self.R = R
        self.r = r
        self.r_s = r_s
        self.B = B
        self.b = b
        m1 = rho*b*(r**2-r_s**2)
        m2 = rho*B*(R**2-r**2)
        self.J = 1/2*m1*(r**2+r_s**2)+1/2*m2*(R**2+r**2)

    def torque(self, BETA1, BETA2, BETA_DDOT, n=200):
        b1, b2 = np.radians(BETA1), np.radians(BETA2)
        BETA = np.linspace(b1, b2, n)

        # gimbal speed profile vs. angle
        w_acc = np.sqrt(2 * BETA_DDOT * (BETA - b1))
        w_dec = np.sqrt(2 * BETA_DDOT * (b2 - BETA))
        Beta_dot = np.minimum.reduce([w_acc, w_dec, np.full(n, self.Beta_dot)])

        # output torque of the scissored pair
        Tau = 2 * self.J * self.Theta_dot * Beta_dot * np.sin(BETA)

        # time: dt = dβ / mean speed over each interval
        dBeta = np.diff(BETA)
        w_avg = 0.5 * (Beta_dot[:-1] + Beta_dot[1:])
        time = np.concatenate([[0], np.cumsum(dBeta / w_avg)])

        return Tau, time, BETA

    def plot_torque(self, Tau, time, BETA, ax=None, label=None, xaxis='angle'):
        if ax is None:
            fig, ax = plt.subplots()

        if xaxis == 'angle':
            ax.plot(np.degrees(BETA), Tau, label=label)
            ax.set_xlabel('Gimbal Angle (deg)')
        else:
            ax.plot(time, Tau, label=label)
            ax.set_xlabel('Time (s)')

        ax.set_ylabel('Torque (Nm)')
        ax.grid(True)
        return ax


CMG1 = gyro(gimbalspeed=5, rotospeed=1000, rho=7850, R=0.04, r=0.03, r_s=0.02, B=0.03, b=0.02)
CMG2 = gyro(gimbalspeed=5, rotospeed=1500, rho=7850, R=0.04, r=0.03, r_s=0.02, B=0.03, b=0.02)

fig, ax = plt.subplots()
Tau, time, BETA = CMG1.torque(BETA1=90, BETA2=180, BETA_DDOT=200)
CMG1.plot_torque(Tau, time, BETA, ax=ax, label="90-180", xaxis='time')
#Tau, time, BETA = CMG1.torque(BETA1=0, BETA2=180, BETA_DDOT=200)
#CMG1.plot_torque(Tau, time, BETA, ax=ax, label="0-180", xaxis='time')

ax.legend()
plt.show()


