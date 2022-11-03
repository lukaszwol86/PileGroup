from PileLoads import PileGroup
import pandas as pd
import time

start_time = time.time()

Loads = pd.read_csv('Loads.csv', delimiter=";")
Loads_ULS = Loads[Loads.Type=="ULS"]
Loads_ULS.to_csv('Loads_ULS.csv',index=False,sep =';')

Loads_SLS = Loads[Loads.Type=="SLS"]
Loads_SLS.to_csv('Loads_SLS.csv',index=False,sep =';')

Loads_QP = Loads[Loads['Type']=='SLS_QP']
Loads_QP.to_csv('Loads_QP.csv',index=False,sep =';')

PilesULS = PileGroup()
PilesULS.pile_parameter(Ep=33*10**9,Gp=12*10**9,Ap=109378*10**-6,Jp=58791*10**(4-4*3),m=0,soil_type="cohesive",kd=4615*1000)
PilesULS.loads('Loads_ULS.csv')
#PilesULS.input_piles('piles.csv')
PilesULS.input_piles('final_piles.csv')
PilesULS.piles_matrix()
PilesULS.forces()


PilesSLS = PileGroup()
PilesSLS.pile_parameter(Ep=33*10**9,Gp=12*10**9,Ap=109378*10**-6,Jp=58791*10**(4-4*3),m=0,soil_type="cohesive",kd=9000*1000)
PilesSLS.loads('Loads_SLS.csv')
#PilesSLS.input_piles('piles.csv')
PilesSLS.input_piles('final_piles.csv')
PilesSLS.piles_matrix()
PilesSLS.forces()


PilesQP = PileGroup()
PilesQP.pile_parameter(Ep=33*10**9,Gp=12*10**9,Ap=109378*10**-6,Jp=58791*10**(4-4*3),m=0,soil_type="cohesive",kd=9000*1000)
PilesQP.loads('Loads_QP.csv')
#PilesQP.input_piles('piles.csv')
PilesQP.input_piles('final_piles.csv')
PilesQP.piles_matrix()
PilesQP.forces()

print("--- %s seconds ---" % (time.time() - start_time))

print(PilesULS.normal_forces())
print(PilesSLS.normal_forces())
print(PilesQP.normal_forces())



print(PilesULS.U[0][1])
print(PilesULS.deformations())
print(PilesSLS.deformations())
print(PilesQP.deformations())