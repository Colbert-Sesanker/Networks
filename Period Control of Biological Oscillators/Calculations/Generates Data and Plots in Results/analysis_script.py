# Generates plots and visulizations
from analysis import *


##################################################################################
#########################       Negative Feedback       ##########################
##################################################################################
# Repressilator
from models import repressilator
R = calc_ens(repressilator, Fs=2, scale=1, steps=300, eps=3e-1)
ens, params = plot_ens(repressilator, Fs=2)
ensemble(*R, scale=5000)
cd = control_direction('repressilator_ens')[0]
control_plot(repressilator, [0, 1], cd[9], var='A')
oscillate_test(repressilator, Fs=2.5, var='A')
delta(params, ens[-1])

# Goodwin
from models import goodwin
G = calc_ens(goodwin, Fs=2.5, scale=1, eps=3e-1)
params = get_parameters(goodwin)
ens, params= plot_ens(goodwin, Fs=2.5, 
                      ens_file='goodwin_dec_ens_final')
ensemble(*G, scale=1,steps=500)
eigen = control_direction('goodwin_ens_inc_final')[0]
params_G = params_array(params)
cd = (ens.copy()[-1] - params_G)
control_plot(goodwin, [0, 1],  cd, var='x')
oscillate_test(goodwin, Fs=2.5, var='x',time=600)
delta(params, ens[100])
# Metabolator


##################################################################################
#########################       Positive Feedback       ##########################
##################################################################################

# Repressilator with auto-activation
from models import repressilator_positive
RP_net = repressilator_positive.Network
RP = calc_ens(repressilator_positive, Fs=2, 
              scale=15, steps=2000, eps=3e-1)
ens, params = plot_ens(repressilator_positive,
                       Fs=2,ens_file='repressilator_positive_dec_ens_final')

ensemble(*RP, scale=5000)
eigen = control_direction('repressilator_positive_ens')[0]
params_RP = params_array(params)
cd = (ens.copy()[-1] - params_RP)
control_plot(repressilator_positive, [0, 1],  cd, var='A')
oscillate_test(repressilator_positive, Fs=2.5, var='A')
delta(params, ens[-1])

# Substrate Depletion
from models import SubstrateDepletion
SD = calc_ens(SubstrateDepletion, scale=5, steps=1000)
ens, params = plot_ens(SubstrateDepletion, 
                       ens_file= 'SubstrateDepletion_inc_ens_final' )
params_SD = params_array(params)
cd = (ens.copy()[-1] - params_SD)
control_plot(SubstrateDepletion, [0, 1],  cd, var='X')
ensemble(*SD)
eigen = control_direction('SubstrateDepletion_inc_ens_final')

# Calcium Oscilations
from models import calcium_oscillations
CO = calc_ens(calcium_oscillations, Fs=2.5, 
              scale=5,steps=300,)
ens, params=plot_ens(calcium_oscillations,
                     ens_file='calcium_oscillations_inc_ens_final', Fs=2.5 )
ensemble(*CO, scale=5)
cd = control_direction('calcium_oscillations_ens')[0]
params_CO = params_array(params)
cd = (ens.copy()[-1] - params_CO)
control_plot(calcium_oscillations, [0, 1], cd, var='x')
oscillate_test(calcium_oscillations, Fs=2.5, var='x')

# Calcium Oscilations 2
from models import calcium_oscillations_2
CO2 = calc_ens(calcium_oscillations_2, Fs=2.5, 
              scale=1,steps=300, eps=3e-1)
ens, params=plot_ens(calcium_oscillations_2, Fs=2.5 )
ensemble(*CO2, scale=5000)
cd = control_direction('calcium_oscillations_2__ens')[0]
control_plot(calcium_oscillations, [0, 1], cd[9], var='x')
oscillate_test(calcium_oscillations_2, Fs=2.5, var='x')
