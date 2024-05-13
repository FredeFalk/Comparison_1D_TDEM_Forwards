import numpy as np
import time
from SimPEG import maps
import SimPEG.electromagnetics.time_domain as tdem
from timeit import default_timer as tictoc

if usestep == 1:
    general_waveform = tdem.sources.StepOffWaveform()
    
else:
        general_waveform = tdem.sources.PiecewiseLinearWaveform(
            times=waveform_times, currents=waveform_current
        )
        
a=time.time()

conductivities = np.array(conductivities)
thicknesses = np.array(thicknesses)
cond_shape = np.shape(conductivities)

NDIM = conductivities.ndim
NM = int(nmodz)

receiver_location = receiver
receiver_orientation = "z"  

receiver_list = [
    tdem.receivers.PointMagneticFluxTimeDerivative(
        receiver_location, times, orientation=receiver_orientation
    )
]

source_location = source
current_amplitude = 1.0

source_list = []

# General Waveform
source_list.append(
    tdem.sources.LineCurrent(
        receiver_list=receiver_list,
        location=source_location,
        current=current_amplitude,
        waveform = general_waveform,
    )
)

survey = tdem.Survey(source_list)

#######################################################################
# Define the Forward Simulation and Predict Data
# ----------------------------------------------
b = time.time()
# Define the simulation

NT = np.size(times)
dpred = np.zeros([NM,NT])
int_time = np.zeros([NM])

simulation = tdem.Simulation1DLayered(
    survey=survey, thicknesses=np.array([]), sigmaMap=maps.IdentityMap(nP=1))

t1 = tictoc()
for i in range(NM):

    if NM > 1:
        n_layer = Nlayers[i]
        
        if n_layer > 1:

            ind1 = 0
            ind2 = int(n_layer-1)
            ind3 = int(n_layer)

            ts = thicknesses[i][ind1:ind2]
            sigma_model = conductivities[i][ind1:ind3]

        else:
            ts = []
            sigma_model = conductivities[i][0]
    else:
        n_layer = int(Nlayers)
        if n_layer > 1:
            ts = thicknesses[0:(n_layer-1)]
            sigma_model = conductivities[0:n_layer]
        else:
            ts = []
            sigma_model = conductivities[0]
    
    # Define a mapping for conductivities
#    model_mapping = maps.IdentityMap(nP=n_layer)
    
    simulation.sigmaMap = maps.IdentityMap(nP=n_layer)
    simulation.thicknesses = ts
    
#    # Thicknesses    
#    simulation = tdem.Simulation1DLayered(
#        survey=survey, thicknesses=ts, sigmaMap=model_mapping,
#    )
        
    # Predict data for a given model
    c=tictoc()
    cdp = simulation.dpred(sigma_model)
    d=tictoc()
    
    dpred[i,:] = cdp
    int_time[i] = d-c


t2 = tictoc()
time_calc_int = np.sum(int_time)
print('Total internal calc time: ',time_calc_int)
time_calc_imp = t2-t1