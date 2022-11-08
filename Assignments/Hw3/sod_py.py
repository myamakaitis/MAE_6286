import numpy as np
import matplotlib.pyplot as pyp


#%%
# -------- Write code --------
def cSound(gamma, Rspecific, TempK):
    return np.sqrt(gamma * Rspecific * TempK)


def flux(u1, u2, u3, gm1 = (1.4-1)):
    f1 = u2
    f2 = u2**2 / u1 + gm1*(u3 - 0.5*(u2**2 / u1))
    f3 = u2/u1 * (u3 + gm1*(u3 - 0.5*(u2**2/u1)))

    return np.array([f1, f2, f3])


class ShockTube:
    def __init__(self, N, Xi, Xe, Xdiaphragm, P0, R0, Ps, Rs, gamma):
        self.X = np.linspace(Xi, Xe, N)
        self.Xdiaphragm = Xdiaphragm
        self.dx = (Xe - Xi) / (N - 1) # I spent like 2 hours bugfixing bc I switched Xi and Xe so my dx was negative
                                      # I was looking all over my flux terms for sign errors and it was right here
        self.U = np.empty((3, N), dtype=np.float64)
        self.F = np.empty((3, N), dtype=np.float64)

        self.U_RM = np.empty((3, (N - 1)), dtype=np.float64)
        self.F_RM = np.empty((3, (N - 1)), dtype=np.float64)

        self.gamma = gamma
        self.gm1 = self.gamma - 1

        self.t = 0

        self.U[0, self.X < Xdiaphragm] = R0
        self.U[0, self.X >= Xdiaphragm] = Rs

        self.U[1] = 0  # initial velocity = 0, therefore momentum = 0

        e0 = P0 / ((self.gm1) * R0)
        es = Ps / ((self.gm1) * Rs)
        self.U[2, self.X < Xdiaphragm] = R0 * e0
        self.U[2, self.X >= Xdiaphragm] = Rs * es

    def BC(self, Flux):
        # Apply fluxes to boundary cells
        self.U[:, 0] -= (dt / self.dx) * (Flux[:, 0] - 0)  # Flux at 1/2 - 0 flux at -1/2
        self.U[:, -1] -= (dt / self.dx) * (0 - Flux[:, -1])  # 0 Flux at N + 1/2 - flux at N - 1/2

        # Velocity = 0 at boundary cells
        # There is a momentum flux which enforces this bc the walls of the tube can exert a pressure
        self.U[1, 0] = 0
        self.U[1, -1] = 0

    def f(self):
        # caclulates flux of each cell
        # self.F[0] = self.U[1]
        # self.F[1] = self.U[1] ** 2 / self.U[0] + self.gm1 * (self.U[2] + .5 * (self.U[1] ** 2 / self.U[0]))
        # self.F[2] = ((self.U[2] + self.gm1 * (self.U[2] - .5 * self.U[1] ** 2 / self.U[0]))
        #              * self.U[1] ** 2 / self.U[0])
        self.F = flux(self.U[0], self.U[1], self.U[2], gm1=self.gm1)

    def f_RM(self):
        # caclulates flux for each 1/2 value in Richtmeyer method
        # This is a huge block of duplicate code that I'm ashamed of
        # My folly is to never want to pass entire arrays to functions
        # self.F_RM[0] = self.U_RM[1]
        # self.F_RM[1] = self.U_RM[1] ** 2 / self.U_RM[0] + self.gm1 * (self.U_RM[2] + .5 * (self.U_RM[1] ** 2 / self.U_RM[0]))
        # self.F_RM[2] = ((self.U_RM[2] + self.gm1 * (self.U_RM[2] - .5 * self.U_RM[1] ** 2 / self.U_RM[0]))
        #                 * self.U_RM[1] ** 2 / self.U_RM[0])
        self.F_RM = flux(self.U_RM[0], self.U_RM[1], self.U_RM[2], gm1=self.gm1)

    def LaxFriedrichsFlux(self, dt):
        self.f()
        # Flux at cell boundaries using lax friedrichs
        self.F_LF = .5 * (self.F[:, :-1] + self.F[:, 1:] - (self.dx / dt) * (self.U[:, 1:] - self.U[:, :-1]))
        # The fact that this ^ this dx/dt is about to be cancelled out annoys me
        # The code is much more readable like this, but it seems wrong

    def GodunovLF(self, dt):
        self.LaxFriedrichsFlux(dt)  # calculates flux using
        self.U[:, 1:-1] -= (dt / self.dx) * (self.F_LF[:, 1:] - self.F_LF[:, :-1])
        self.BC(self.F_LF)  # apply boundary conditions using LF flux

    def Richtmeyer(self, dt):
        self.f()  # calculate flux
        self.U_RM = .5 * (self.U[:, 1:] + self.U[:, :-1]) - (dt / (2 * self.dx)) * (self.F[:, 1:] - self.F[:, :-1])
        # The u now has the values for between the cells
        # self.BC_RM()
        self.f_RM()

        self.U[:, 1:-1] -= (dt / self.dx) * (self.F_RM[:, 1:] - self.F_RM[:, :-1])
        self.BC(self.F_RM)

    def Solve_LF(self, dt, tlim):
        steps = int(np.ceil(tlim / dt))
        for n in range(steps):
            self.GodunovLF(dt)
            self.t += dt

    def Solve_RM(self, dt, tlim):
        steps = int(np.ceil(tlim / dt))
        for n in range(steps):
            self.Richtmeyer(dt)
            self.t += dt

    def Primitives(self, R=287):
        u = self.U[1] / self.U[0]
        rho = self.U[0]
        P = self.gm1 *(self.U[2] - .5*self.U[1]**2/self.U[0])
        T = P / (self.U[0] * R)

        return u, rho, P, T

    def Mach(self, Rgas=287.0):
        u, _, _, T = self.Primitives()
        c = cSound(self.gamma, Rgas, T)
        M = u / c
        return M


dt = 0.0002
Tend = 0.01
Gamma = 1.4

nx = 81

startL = -10
endL = 10
diaphragmL = 0

rL = 1.0 # kg / m^3 Density on left
PL = 100000 # kPA Pressure on left

rR = 0.125 # kg / m^3 Density on right
PR = 10000 # kPA Pressure on right

Tube_LF = ShockTube(nx, startL, endL, diaphragmL, PL, rL, PR, rR, Gamma)
Tube_RM = ShockTube(nx, startL, endL, diaphragmL, PL, rL, PR, rR, Gamma)

Tube_LF.Solve_LF(dt, Tend)
Tube_RM.Solve_RM(dt, Tend)

#%%
# fig, axes = pyp.subplots(3, figsize = (7,7), tight_layout = True)
# for n in range(100):
#     Tube_LF.GodunovLF(dt)
# for n in range(25):
#     Tube_LF.GodunovLF(dt)
#     if n%5 == 0:
#         for i, ax in enumerate(axes):
#             ax.set_title(f"u{i}")
#             ax.plot(Tube_LF.X, Tube_LF.U[i], label=n+100)
#             ax.legend()
# fig.show()


u_LF, rho_LF, P_LF, T_LF = Tube_LF.Primitives()
u_RM, rho_RM, P_RM, T_RM = Tube_RM.Primitives()


# -------- Write code --------
from sod import analytical_solution
X = np.linspace(startL, endL, nx)
rho_anal, u_anal, P_anal = analytical_solution(Tend, X, (rL, 0, PL), (rR, 0, PR))

fig, (axr, axu, axP) = pyp.subplots(3, figsize = (7,7), tight_layout = True)


axr.plot(X, rho_anal, color='k', linestyle =':', label='Analytical')
axr.plot(Tube_LF.X, rho_LF, color='darkorange', label="Lax Friedrichs")
axr.plot(Tube_RM.X, rho_RM, color='rebeccapurple', label="Richtmeyer")
axr.set_ylabel("Density [kg / m^3]")
axr.legend()


axu.plot(X, u_anal, color='k', linestyle =':')
axu.plot(Tube_LF.X, u_LF, color='darkorange')
axu.plot(Tube_RM.X, u_RM, color='rebeccapurple')
axu.set_ylabel("Velocity [m/s]")


axP.plot(X, P_anal/1e3, color='k', linestyle =':')
axP.plot(Tube_LF.X, P_LF/1e3, color='darkorange')
axP.plot(Tube_RM.X, P_RM/1e3, color='rebeccapurple')
axP.set_ylabel("Pressure [kPa]")

axP.set_xlabel("X [m]")

fig.show()
