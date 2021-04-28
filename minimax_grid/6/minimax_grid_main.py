from modi_leastsq import main_freq
from tau_grid import main_time
from time2freq_transform_grid import main_trans

Npoint=6
R=90
if __name__=='__main__':
    print("Run_freq_grid")
    main_freq(Npoint,R)
    print("Run_time_grid")
    main_time(Npoint,R)
    print("Run_transform_grid")
    main_trans(R)
    print("Finish_minimax_grid")
