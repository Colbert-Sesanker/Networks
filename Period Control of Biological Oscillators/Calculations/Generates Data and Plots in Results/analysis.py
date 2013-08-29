# Provides support for the analysis script

import os
import scipy
import scipy.io
from ReactionNetworks import *
from models import *
from pylab  import *
from oct2py import octave as oc

def ensemble(m, params, ens_file, hess, scale=1.0, steps=200, save_hours = 1/(60*3.0)): 
    Ensembles.ensemble(m, asarray(params), hess=hess, steps = steps, 
                       save_to=ens_file,  periodic=True, step_scale=scale, 
                       save_hours=save_hours)

def get_parameters(model):
    reload(model)
    net = model.network
    return(net.GetParameters())

def calc_ens(model, Fs=10, relativeScale=True, scale=1.0, 
             steps=200, eps=1e-1, time=None):
  reload(model)
  net             = model.network
  if time is None:
      time_interval   = model.int_time
      time            = time_interval[1] - time_interval[0] 
  else:
      time_interval   = [0,time]    
  num_data_points = time*Fs
  params          = net.GetParameters()
  model_name      = net.get_id()
  ens_file        = model_name + '_ens'
  
  # WARNING: Network name must match experiment name
  fake_expt = Experiment('synthetic_expt')
  def create_expt(data, net): fake_expt.SetData({model_name: data})

  perfect_data    = PerfectData.discrete_data(net, params,   
                                              num_data_points, time_interval)
  create_expt(perfect_data, net)
  m = Model([fake_expt], [net])  
 
  hess = m.oscillate_hessian(params, eps, relativeScale=relativeScale)  
  try:
      ensemble(m, params, ens_file, hess, scale, steps)
  except:
      print "trajectory blew up:" 
      return m, params, ens_file, hess
  return m, params, ens_file, hess

def plot_ens(model, ens_file=None, Fs=10):
  reload(model)
  net             = model.network
  time_interval   = model.int_time
  time            = time_interval[1] - time_interval[0]  
  num_data_points = time*Fs
  variables       = net.dynamicVars.keys()
  params          = net.GetParameters()  
  if ens_file is None:  
      model_name      = net.get_id()
      ens_file    = model_name + '_ens'

  import ReactionNetworks.Utility as Utility
  ens =  Utility.load(ens_file)[0]  
  prn = int(ceil(len(ens) / 7))                        
  pruned_ens = asarray(ens[::prn])
  times = linspace(0, time, num_data_points)
  freqs = linspace(0, Fs/2, ceil(num_data_points / 2))
  ens_trajs = Ensembles.ensemble_trajs(net, times, pruned_ens)

  def get_var_traj(var):
      return([traj.get_var_traj(var) for traj in ens_trajs])

  def get_var_spec(var):
      return([oc.pow_spec(traj.get_var_traj(var))[0] for traj in ens_trajs])

  ion(); draw();
  def plot_traj(times, var_traj_list,  inc):    
    for traj in range(0, len(var_traj_list), inc):        
          plot(times, var_traj_list[traj])
  
  for var in variables:
      var_traj = get_var_traj(var)
      figure()
      suptitle(var + ' ensemble (Abundance)')      
      plot_traj(times, var_traj[:int(100*Fs)],  1)
      xlabel('time (seconds)')
      ylabel(var)

      figure()
      suptitle(var + ' ensemble (Normalized Power Spectrum)')    
      spec = get_var_spec(var)        
      plot_traj(freqs, spec,  1)      
      xlabel('Frequency (Hz)')
      ylabel('Power')
  return asarray(ens), params
############# Frequency Control ###############

def control_direction(ens_file):
    import ReactionNetworks.Utility as Utility
    ens =  Utility.load(ens_file)[0]
    eigs_PCA = Ensembles.PCA_eig(ens)
    eig_vecs_PCA = eigs_PCA[1]
    return eig_vecs_PCA, eigs_PCA 

def params_array(params):
    return([param[1] for param in params.items()])

def delta(params, centroid):
    keys = params.keys()
    delta = ((asarray(centroid)- asarray(params)) /
              asarray(params))
    return(zip(keys, delta))

    # Slide parameters along control parameter
    # net_params is a keyed list
def get_controled_net_params(net_params, knob, control_parameter):
    controled_net_params  = net_params.deepcopy()
    control_parameter = array(control_parameter)*knob
    for cp_indx, cp in enumerate(control_parameter):
        controled_net_params[cp_indx] += cp
    return controled_net_params

def control_plot(model, knob_range, control_parameter, 
                 var,samples=5, Fs=10):
    reload(model)
    net             = model.network.copy()
    time_interval   = model.int_time
    time            = time_interval[1] - time_interval[0]  
    num_data_points = time*Fs
    params          = net.GetParameters().deepcopy()
    times = linspace(0, time, num_data_points)
    freqs = linspace(0, Fs/2, ceil(num_data_points / 2))
    cp    = control_parameter.copy()
    
    knobs = linspace(knob_range[0], knob_range[1], samples)
    for knob in knobs:            
        cp    = get_controled_net_params(params, knob, 
                                         control_parameter)
        trajs = Dynamics.integrate(net, times,
                                   params=cp.deepcopy(), 
                                   fill_traj=False)         
        ion(); draw;        
        figure(1)
        suptitle(var + ' Abundance')
        xlabel('time (seconds)')
        ylabel(var)
        plot(times, trajs.get_var_traj(var))
        figure(2)
        suptitle(var + ' Normalized Power Spectrum')
        xlabel('Frequency (Hz)')
        ylabel('Power')
        plot(freqs, oc.pow_spec(trajs.get_var_traj(var))[0])        
    
# Quick Test to see if network is oscillating
def oscillate_test(model, Fs, var, time=None):
    reload(model)
    net             = model.network.copy()    
    if time is None:
        time_interval   = model.int_time
        time            = time_interval[1] - time_interval[0]  
    num_data_points = time*Fs
    params          = net.GetParameters().deepcopy()
    times = linspace(0, time, num_data_points)
    freqs = linspace(0, Fs/2, ceil(num_data_points / 2))
    trajs = Dynamics.integrate(net, times,
                               params=params, 
                               fill_traj=False)
    ion(); draw;        
    figure(1)
    suptitle(var + ' Abundance')
    xlabel('time (seconds)')
    ylabel(var)
    plot(times, trajs.get_var_traj(var))
    figure(2)
    suptitle(var + ' Normalized Power Spectrum')
    xlabel('Frequency (Hz)')
    ylabel('Power')
    plot(freqs, oc.pow_spec(trajs.get_var_traj(var))[0]) 

#control_plot(params,[1,100], eig_vecs_PCA[9])

#get_controled_net_params(params, 9, 
# eig_vecs_PCA[0])

















