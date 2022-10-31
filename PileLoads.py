import numpy as np
import pandas as pd
import math


class PileGroup():
    np.set_printoptions(formatter={'float': '{: 0.2f}'.format})
    def pile_parameter(self,Ep,Gp,Ap,Jp,m,soil_type,**kwarg):
        self.Ep = Ep
        self.Gp = Gp
        self.Ap = Ap
        self.Jp = Jp
        self.m = m
        self.soil_type = soil_type
        if soil_type=="cohesive":
            self.kd= kwarg["kd"]
        elif soil_type=="friction":
            self.nh = kwarg["nh"]


    def loads(self,input_loads):
        self.Loads = pd.read_csv(input_loads, delimiter=";")
        self.R = []
        for index, row in self.Loads.iterrows():
            r = np.array([row.Fx * 1000, row.Fy * 1000, row.Fz * 1000, row.Mx * 1000, row.My * 1000, row.Mz * 1000])
            self.R.append(r)

    def single_load(self,input_loads):
        self.Loads = pd.read_csv(input_loads, delimiter=";")
        self.R = []
        row = self.Loads.loc[1]
        r = np.array([row.Fx * 1000, row.Fy * 1000, row.Fz * 1000, row.Mx * 1000, row.My * 1000, row.Mz * 1000])
        self.R.append(r)

    def input_piles(self,input_piles_pos):
        self.piles = pd.read_csv(input_piles_pos, delimiter=";")

    def import_piles(self,piles_df):
        self.piles = piles_df[piles_df['Acti']==1]


    def piles_matrix(self):
        self.Kp = []
        for index, row in self.piles.iterrows():
            L = row.L
            k = np.zeros((6, 6))
            if self.soil_type == "no":
                k[0, 0] = self.m * 3 * self.Ep * self.Jp / L ** 3
                k[1, 1] = k[0, 0]
                k[2, 2] = self.Ap * self.Ep / L
                k[3, 3] = self.m * 3 * self.Ep * self.Jp / L
                k[4, 4] = k[3, 3]
                k[5, 5] = self.m * self.Gp * self.Jp / L
            elif self.soil_type == "cohesive":
                Le = (4 * self.Ep * self.Jp / (self.kd*row.redu)) ** (1 / 4)
                k[0, 0] = (self.m + 1) * 2 * self.Ep * self.Jp / Le ** 3
                k[1, 1] = k[0, 0]
                k[2, 2] = self.Ap * self.Ep / L
                k[3, 3] = self.m * 2 * self.Ep * self.Jp / Le
                k[4, 4] = k[3, 3]
                k[5, 5] = self.m * self.Gp * self.Jp / L
                k[0, 4] = self.m * 2 * self.Ep * self.Jp / Le ** 2
                k[4, 0] = k[0, 4]
                k[1, 3] = -k[0, 4]
                k[3, 1] = -k[0, 4]
                # print(f"Pile {index+1} Le= {Le: 0.2f}")

            elif self.soil_type == "friction":
                Li = 1.8 * (self.Ep * self.Jp / (self.nh*row.redu)) ** (1 / 5)
                k[0, 0] = (3 * self.m + 1) * 3 * self.Ep * self.Jp / Li ** 3
                k[1, 1] = k[0, 0]
                k[2, 2] = self.Ap * self.Ep / L
                k[3, 3] = self.m * 4 * self.Ep * self.Jp / Li
                k[4, 4] = k[3, 3]
                k[5, 5] = self.m * self.Gp * self.Jp / L
                k[0, 4] = self.m * 6 * self.Ep * self.Jp / Li ** 2
                k[4, 0] = k[0, 4]
                k[1, 3] = -k[0, 4]
                k[3, 1] = -k[0, 4]

            self.Kp.append(k)

        self.Dp = []
        self.Cp = []
        self.Api = []
        self.S = np.zeros((6, 6))
        i=0
        for index, row in self.piles.iterrows():
            cp = np.array([[1, 0, 0, 0, 0, 0],
                           [0, 1, 0, 0, 0, 0],
                           [0, 0, 1, 0, 0, 0],
                           [0, -row.z, row.y, self.m, 0, 0],
                           [row.z, 0, -row.x, 0, self.m, 0],
                           [-row.y, row.x, 0, 0, 0, self.m]])
            beta = math.atan(1 / row.N)
            alfa = math.radians(row.Angle)

            ap = np.array([[math.cos(beta) * math.cos(alfa), -math.sin(alfa), math.sin(beta) * math.cos(alfa), 0, 0, 0],
                           [math.cos(beta) * math.sin(alfa), math.cos(alfa), math.sin(beta) * math.sin(alfa), 0, 0, 0],
                           [-math.sin(beta), 0, math.cos(beta), 0, 0, 0],
                           [0, 0, 0, math.cos(beta) * math.cos(alfa), -math.sin(alfa), math.sin(beta) * math.cos(alfa)],
                           [0, 0, 0, math.cos(beta) * math.sin(alfa), math.cos(alfa), math.sin(beta) * math.sin(alfa)],
                           [0, 0, 0, -math.sin(beta), 0, math.cos(beta)]])

            dp = np.matmul(cp, ap)
            sp = np.matmul(np.matmul(dp, self.Kp[i]), dp.T)

            self.Dp.append(dp)
            self.S = np.add(self.S, sp)

            self.Api.append(ap)
            self.Cp.append(cp)
            i+=1

    def forces(self):
        self.U = []
        self.DispP = np.zeros((len(self.Loads), len(self.piles), 6))
        self.ForcsP = np.zeros((len(self.Loads), len(self.piles), 6))
        for i, reaction in enumerate(self.R):
            u = np.matmul(np.linalg.inv(self.S), reaction)
            self.U.append(u)
            for n, dp in enumerate(self.Dp):
                disp = np.matmul(dp.T, u)
                force = np.matmul(self.Kp[n], disp) / 1000
                # print(disp)
                self.DispP[i, n] = disp
                self.ForcsP[i, n] = force

    def single_normal_forces(self):
        N_P_LC = np.zeros(len(self.piles))

        for pile in range(len(self.piles)):
                N_P_LC[ pile] = self.ForcsP[0, pile, 2]

        return N_P_LC

    def single_deformation(self):
        deform = ((self.U[0][0]) ** 2 + (self.U[0][1]) ** 2) ** 0.5 * 1000
        return deform


    def normal_forces(self):
        N_P_LC = np.zeros((len(self.Loads), len(self.piles)))

        for pile in range(len(self.piles)):
            for lc in range(len(self.Loads)):
                N_P_LC[lc, pile] = self.ForcsP[lc, pile, 2]

        return N_P_LC

    def deformations(self):
        deform = np.zeros((len(self.Loads), 4))

        for lc in range(len(self.Loads)):
            deform[lc, 0] = self.U[lc][0]*1000
            deform[lc, 1] = self.U[lc][1]*1000
            deform[lc, 2] = self.U[lc][2]*1000
            deform[lc, 3] = ((self.U[lc][0])**2+(self.U[lc][1])**2)**0.5*1000
        return deform


if __name__ == "__main__":
    Piles = PileGroup()
    Piles.pile_parameter(Ep=33*10**9,Gp=12*10**9,Ap=109378*10**-6,Jp=58791*10**(4-4*3),m=0,soil_type="cohesive",kd=4615*1000)
    Piles.loads('Loads.csv')
    Piles.input_piles('piles.csv')
    Piles.piles_matrix()
    Piles.forces()



    print(Piles.normal_forces())








