
'''
Script to access simuated GEDI data
'''

##################################
import numpy as np
import h5py
from sys import exit
import matplotlib.pyplot as plt
if __name__ == '__main__':
  import argparse

###################################

class gediData(object):
  '''
  Simulated GEDI data handler
  '''

  def __init__(self,filename=None,minX=-100000000,maxX=100000000,minY=-1000000000,maxY=100000000,outName="teast.h5"):
    '''
    Class initialiser. Calls a function
    to read waveforms between bounds
    and writes to a new file
    '''

    self.doneAnc=False
    self.subsetGEDI(filename,minX,maxX,minY,maxY,outName)


  ###########################################

  def subsetGEDI(self,filename,minX,maxX,minY,maxY,outNamen):
    '''
    Read real GEDI data from file
    '''
    # open file for reading
    f=h5py.File(filename,'r')
    self.beamList=['BEAM0000', 'BEAM0001', 'BEAM0010', 'BEAM0011', 'BEAM0101', 'BEAM0110', 'BEAM1000', 'BEAM1011']
    self.nWaves=0
    fileFormat=None

    # open output file
    outFile=h5py.File(outNamen,'w')

    # set directory list
    self.setRealList()


    # loop over beams
    for b in self.beamList:
      if((b in list(f))==False): # does this exist?
        continue                 # if not, skip it
      elif(('geolocation' in list(f[b]))==False):  # no data in bea,
        continue

      print(b)

      # determine what file version
      if('rxwaveform' in list(f[b])):  # L1B file
        self.subsetL1B(f,outFile,b,minX,maxX,minY,maxY)
        fileFormat="L1B"
      elif ('rh' in list(f[b])):             # L2A file
        self.subsetL2A(f,outFile,b,minX,maxX,minY,maxY)
        fileFormat="L2A"
      elif ('pai' in list(f[b])):             # L2B file  
        self.subsetL2B(f,outFile,b,minX,maxX,minY,maxY)
        fileFormat="L2B"
      elif ('agbd' in list(f[b])):             # L4A file 
        self.subsetL4A(f,outFile,b,minX,maxX,minY,maxY)
        fileFormat="L4A"
      else:
        print('File format not known')
        exit


    # write metadata group if needed
    if((fileFormat=="L2A")|(fileFormat=="L2B")):
      self.metadataL2A(f,outFile)

    f.close()
    outFile.close()
    print("Written to",outNamen)
    return


  ###########################################

  def subsetL4A(self,f,outFile,b,minX,maxX,minY,maxY):
    '''Subset an L4A beam'''

    print('For now this script does not include the model_data information in the subset file')

    # read the coords and determine output
    allLat=np.array(f[b]['lat_lowestmode'])
    allLon=np.array(f[b]['lon_lowestmode'])
    useInd=np.where((allLat>=minY)&(allLat<=maxY)&(allLon>=minX)&(allLon<=maxX))

    if(len(useInd[0])>0):
      useInd=useInd[0]
    else:      # none in here
      return

    # create the beam group
    outFile.create_group(b)

    # loop over all arrays per shot
    for d in self.shotArrListL4A:
      # read array
      jimlad=np.array(f[b][d])[useInd]
      # write subset to a new file
      outFile[b].create_dataset(d,data=jimlad,compression='gzip')

    # string fixed array lengths
    for d in self.strArrListL4A:
      jimlad=np.array(f[b][d])
      # find maximum length
      maxLen=0
      for n in jimlad:
        if(len(n)>maxLen):
          maxLen=len(n)
      asciiList = [n.encode("ascii", "ignore") for n in jimlad]
      outFile[b].create_dataset(d, (len(asciiList),1),'S10', asciiList)

    # geolocation data
    g='geolocation'
    outFile[b].create_group(g)
    for d in self.geoArrListL4A:
      if(d!='surface_type'):
        jimlad=np.array(f[b][g][d])[useInd]
        outFile[b][g].create_dataset(d,data=jimlad,compression='gzip')
      else:
        jimlad=np.array(f[b][g][d])[:,useInd]
        outFile[b][g].create_dataset(d,data=jimlad,compression='gzip')

    if(self.doneAnc==False):
      # ancillary data   
      g='ANCILLARY'
      outFile.create_group(g)
      for d in self.ancArrListL4A:
        jimlad=np.array(f[g][d])
        outFile[g].create_dataset(d,data=jimlad,compression='gzip')
      self.doneAnc=True

    # agbd bits
    g='agbd_prediction'
    outFile[b].create_group(g)
    for d in self.agbdArrListL4A:
      jimlad=np.array(f[b][g][d])
      outFile[b][g].create_dataset(d,data=jimlad,compression='gzip')


    return


  ###########################################

  def subsetL2B(self,f,outFile,b,minX,maxX,minY,maxY):
    '''Subset an L2B beam'''

    # read the coords and determine output
    allLat=np.array(f[b]['geolocation']['lat_lowestmode'])
    allLon=np.array(f[b]['geolocation']['lon_lowestmode'])
    useInd=np.where((allLat>=minY)&(allLat<=maxY)&(allLon>=minX)&(allLon<=maxX))

    if(len(useInd[0])>0):
      useInd=useInd[0]
    else:      # none in here
      return

    # create the beam group
    outFile.create_group(b)

    # loop over all arrays per shot
    for d in self.shotArrListL2B:
      if((d=='rx_sample_start_index')|(d=='tx_sample_start_index')):  # skip these for now
        continue
      # read array
      jimlad=np.array(f[b][d])[useInd]
      # write subset to a new file
      outFile[b].create_dataset(d,data=jimlad,compression='gzip')

    # geolocation
    g='geolocation'
    outFile[b].create_group(g)
    for d in self.geoArrListL2B:
      if(d!='surface_type'):
        jimlad=np.array(f[b][g][d])[useInd]
        outFile[b][g].create_dataset(d,data=jimlad,compression='gzip')
      else:
        jimlad=np.array(f[b][g][d])[:,useInd]
        outFile[b][g].create_dataset(d,data=jimlad,compression='gzip')

    # ancillary data   
    g='ancillary'
    outFile[b].create_group(g)
    for d in self.ancArrListL2B:
      jimlad=np.array(f[b][g][d])
      outFile[b][g].create_dataset(d,data=jimlad,compression='gzip')

    # string fixed array lengths
    for d in self.strArrListL2B:
      jimlad=np.array(f[b][d])
      # find maximum length
      maxLen=0
      for n in jimlad:
        if(len(n)>maxLen):
          maxLen=len(n)
      asciiList = [n.encode("ascii", "ignore") for n in jimlad]
      outFile[b].create_dataset(d, (len(asciiList),1),'S10', asciiList)

    # Variable length array
    self.subsetWaves('pgap_theta_z','rx_sample_start_index','rx_sample_count',useInd,outFile[b],f[b])

    return


  ###########################################

  def subsetL2A(self,f,outFile,b,minX,maxX,minY,maxY):
    '''Subset an L2A beam'''

    # read the coords and determine output
    allLat=np.array(f[b]['lat_lowestmode'])
    allLon=np.array(f[b]['lon_lowestmode'])
    useInd=np.where((allLat>=minY)&(allLat<=maxY)&(allLon>=minX)&(allLon<=maxX))
    
    if(len(useInd[0])>0):
      useInd=useInd[0]
    else:      # none in here
      return

    # create the beam group
    outFile.create_group(b)

    # loop over all arrays per shot
    for d in self.shotArrListL2A:
      # read array
      jimlad=np.array(f[b][d])[useInd]
      # write subset to a new file
      outFile[b].create_dataset(d,data=jimlad,compression='gzip')

    # string fixed array lengths
    for d in self.strArrListL2A:
      jimlad=np.array(f[b][d])
      # find maximum length
      maxLen=0
      for n in jimlad:
        if(len(n)>maxLen):
          maxLen=len(n)
      asciiList = [n.encode("ascii", "ignore") for n in jimlad]
      outFile[b].create_dataset(d, (len(asciiList),1),'S10', asciiList)

    # geolocation data
    g='geolocation'
    outFile[b].create_group(g)
    for d in self.geoArrListL2A:
      if(d!='surface_type'):
        jimlad=np.array(f[b][g][d])[useInd]
        outFile[b][g].create_dataset(d,data=jimlad,compression='gzip')
      else:
        jimlad=np.array(f[b][g][d])[:,useInd]
        outFile[b][g].create_dataset(d,data=jimlad,compression='gzip')

    # ancillary data   
    g='ancillary'      
    outFile[b].create_group(g)
    for d in self.ancArrListL2A:
      jimlad=np.array(f[b][g][d])
      outFile[b][g].create_dataset(d,data=jimlad,compression='gzip')

    return


  ###########################################

  def subsetL1B(self,f,outFile,b,minX,maxX,minY,maxY):
    '''Subset an L1B beam'''

    # read the coords and determine output
    allLat=(np.array(f[b]['geolocation']['latitude_bin0'])+np.array(f[b]['geolocation']['latitude_lastbin']))/2.0
    allLon=(np.array(f[b]['geolocation']['longitude_bin0'])+np.array(f[b]['geolocation']['longitude_lastbin']))/2.0
    useInd=np.where((allLat>=minY)&(allLat<=maxY)&(allLon>=minX)&(allLon<=maxX))

    if(len(useInd[0])>0):
      useInd=useInd[0]
    else:      # none in here
      return

    # create the beam group
    outFile.create_group(b)

    # loop over all arrays per shot
    for d in self.shotArrListL1B:
      if((d=='rx_sample_start_index')|(d=='tx_sample_start_index')):  # skip these for now
        continue
      # read array
      jimlad=np.array(f[b][d])[useInd]
      # write subset to a new file
      outFile[b].create_dataset(d,data=jimlad,compression='gzip')

    # for txwaveform and rxwaveform, we must read start/stop indices
    self.subsetWaves('rxwaveform','rx_sample_start_index','rx_sample_count',useInd,outFile[b],f[b])
    self.subsetWaves('txwaveform','tx_sample_start_index','tx_sample_count',useInd,outFile[b],f[b])

    # geolocation data
    g='geolocation'
    outFile[b].create_group(g)
    for d in self.geoArrListL1B:
      if(d!='surface_type'):
        jimlad=np.array(f[b][g][d])[useInd]
        outFile[b][g].create_dataset(d,data=jimlad,compression='gzip')
      else:
        jimlad=np.array(f[b][g][d])[:,useInd]
        outFile[b][g].create_dataset(d,data=jimlad,compression='gzip')

    # ancillary data
    g='ancillary'
    outFile[b].create_group(g)
    for d in self.ancArrListL1B:
      jimlad=np.array(f[b][g][d])
      outFile[b][g].create_dataset(d,data=jimlad,compression='gzip')

    # geophys_corr
    g='geophys_corr'
    outFile[b].create_group(g)
    for d in self.corArrListL1B:
      jimlad=np.array(f[b][g][d])[useInd]
      outFile[b][g].create_dataset(d,data=jimlad,compression='gzip')

    return


  ###########################################

  def subsetWaves(self,waveName,indName,lenName,useInd,outFile,f):
    '''Subset and write the RX or TX waveforms'''

    # read indices
    startInds=np.array(f[indName])[useInd]
    lenInds=np.array(f[lenName])[useInd]
    totBins=np.sum(lenInds)
    waveform=np.empty(totBins,dtype=np.float32)
    newInds=np.zeros(len(useInd),dtype=np.uint64)

    # read raw data and repack
    jimlad=np.array(f[waveName])
    lastInd=0
    for i in range(0,startInds.shape[0]):
      newInds[i]=lastInd
      waveform[lastInd:lastInd+lenInds[i]]=jimlad[startInds[i]:startInds[i]+lenInds[i]]
      lastInd=lastInd+lenInds[i]

    # write data
    outFile.create_dataset(waveName,data=waveform,compression='gzip')
    outFile.create_dataset(indName,data=newInds,compression='gzip')

    return


  ###########################################

  def metadataL2A(self,f,outFile):
    '''Write metadata group for L2A file'''

    g='METADATA'
    outFile.create_group(g)
    outFile[g].create_group('DatasetIdentification')

    return

  ###########################################

  def setRealList(self):
    '''Set list of all data within a real GEDI file'''

    # L1B files
    # arrays with one element per shot
    self.shotArrListL1B=['all_samples_sum', 'beam', 'channel', 'delta_time', 'master_frac', 'master_int',\
                        'noise_mean_corrected', 'noise_stddev_corrected', 'nsemean_even', 'nsemean_odd',\
                        'rx_energy', 'rx_offset', 'rx_open', 'rx_sample_count',\
                        'selection_stretchers_x', 'selection_stretchers_y', 'shot_number', 'stale_return_flag',\
                        'th_left_used', 'tx_egamplitude', 'tx_egamplitude_error', 'tx_egbias', 'tx_egbias_error',\
                        'tx_egflag', 'tx_eggamma', 'tx_eggamma_error', 'tx_egsigma', 'tx_egsigma_error', 'tx_gloc',\
                        'tx_gloc_error', 'tx_pulseflag', 'tx_sample_count',]

    self.geoArrListL1B=['altitude_instrument', 'altitude_instrument_error', 'bounce_time_offset_bin0', 'bounce_time_offset_bin0_error',\
                       'bounce_time_offset_lastbin', 'bounce_time_offset_lastbin_error', 'degrade', 'delta_time',\
                       'digital_elevation_model', 'elevation_bin0', 'elevation_bin0_error', 'elevation_lastbin',\
                       'elevation_lastbin_error', 'latitude_bin0', 'latitude_bin0_error', 'latitude_instrument',\
                       'latitude_instrument_error', 'latitude_lastbin', 'latitude_lastbin_error', 'local_beam_azimuth',\
                       'local_beam_azimuth_error', 'local_beam_elevation', 'local_beam_elevation_error', 'longitude_bin0',\
                       'longitude_bin0_error', 'longitude_instrument', 'longitude_instrument_error', 'longitude_lastbin',\
                       'longitude_lastbin_error', 'mean_sea_surface', 'neutat_delay_derivative_bin0', 'neutat_delay_derivative_lastbin',\
                       'neutat_delay_total_bin0', 'neutat_delay_total_lastbin', 'range_bias_correction', 'shot_number',\
                       'solar_azimuth', 'solar_elevation', 'surface_type']

    self.ancArrListL1B=['master_time_epoch', 'mean_samples', 'smoothing_width']

    self.corArrListL1B=['delta_time', 'dynamic_atmosphere_correction', 'geoid', 'tide_earth', 'tide_load', 'tide_ocean',\
                       'tide_ocean_pole', 'tide_pole']


    ## L2A files
    # element per shot
    self.shotArrListL2A=['rh','beam','channel','degrade_flag','delta_time','digital_elevation_model','digital_elevation_model_srtm',\
                       'elev_highestreturn','elev_lowestmode','elevation_bias_flag','elevation_bin0_error','energy_total',\
                       'lat_highestreturn', 'lat_lowestmode', 'latitude_bin0_error', 'lon_highestreturn', 'lon_lowestmode',\
                       'longitude_bin0_error','master_frac','master_int','mean_sea_surface','num_detectedmodes','quality_flag',\
                       'selected_algorithm','selected_mode','selected_mode_flag','sensitivity','shot_number','solar_azimuth',\
                       'solar_elevation','surface_flag']

    # geolocation. Element per shot
    self.geoArrListL2A=['elev_highestreturn_a1', 'elev_highestreturn_a2', 'elev_highestreturn_a3', 'elev_highestreturn_a4', \
                       'elev_highestreturn_a5', 'elev_highestreturn_a6', 'elev_lowestmode_a1', 'elev_lowestmode_a2', 'elev_lowestmode_a3',\
                       'elev_lowestmode_a4', 'elev_lowestmode_a5', 'elev_lowestmode_a6', 'elev_lowestreturn_a1', 'elev_lowestreturn_a2',\
                       'elev_lowestreturn_a3', 'elev_lowestreturn_a4', 'elev_lowestreturn_a5', 'elev_lowestreturn_a6', 'elevation_1gfit',\
                       'elevs_allmodes_a1', 'elevs_allmodes_a2', 'elevs_allmodes_a3', 'elevs_allmodes_a4', 'elevs_allmodes_a5', \
                       'elevs_allmodes_a6', 'energy_lowestmode_a1', 'energy_lowestmode_a2', 'energy_lowestmode_a3', 'energy_lowestmode_a4',\
                       'energy_lowestmode_a5', 'energy_lowestmode_a6', 'lat_highestreturn_a1', 'lat_highestreturn_a2', 'lat_highestreturn_a3', \
                       'lat_highestreturn_a4', 'lat_highestreturn_a5', 'lat_highestreturn_a6', 'lat_lowestmode_a1', 'lat_lowestmode_a2',\
                       'lat_lowestmode_a3', 'lat_lowestmode_a4', 'lat_lowestmode_a5', 'lat_lowestmode_a6', 'lat_lowestreturn_a1', \
                       'lat_lowestreturn_a2', 'lat_lowestreturn_a3', 'lat_lowestreturn_a4', 'lat_lowestreturn_a5', 'lat_lowestreturn_a6',\
                       'latitude_1gfit', 'lats_allmodes_a1', 'lats_allmodes_a2', 'lats_allmodes_a3', 'lats_allmodes_a4', 'lats_allmodes_a5',\
                       'lats_allmodes_a6', 'lon_highestreturn_a1', 'lon_highestreturn_a2', 'lon_highestreturn_a3', 'lon_highestreturn_a4',\
                       'lon_highestreturn_a5', 'lon_highestreturn_a6', 'lon_lowestmode_a1', 'lon_lowestmode_a2', 'lon_lowestmode_a3', \
                       'lon_lowestmode_a4', 'lon_lowestmode_a5', 'lon_lowestmode_a6', 'lon_lowestreturn_a1', 'lon_lowestreturn_a2', \
                       'lon_lowestreturn_a3', 'lon_lowestreturn_a4', 'lon_lowestreturn_a5', 'lon_lowestreturn_a6', 'longitude_1gfit',\
                       'lons_allmodes_a1', 'lons_allmodes_a2', 'lons_allmodes_a3', 'lons_allmodes_a4', 'lons_allmodes_a5', 'lons_allmodes_a6',\
                       'num_detectedmodes_a1', 'num_detectedmodes_a2', 'num_detectedmodes_a3', 'num_detectedmodes_a4', 'num_detectedmodes_a5',\
                       'num_detectedmodes_a6', 'quality_flag_a1', 'quality_flag_a2', 'quality_flag_a3', 'quality_flag_a4', 'quality_flag_a5',\
                       'quality_flag_a6', 'rh_a1', 'rh_a2', 'rh_a3', 'rh_a4', 'rh_a5', 'rh_a6', 'sensitivity_a1', 'sensitivity_a2', \
                       'sensitivity_a3', 'sensitivity_a4', 'sensitivity_a5', 'sensitivity_a6', 'shot_number', 'stale_return_flag']

    # Fixed number of elements
    self.strArrListL2A=['land_cover_data','rx_1gaussfit','rx_assess','rx_processing_a1','rx_processing_a2','rx_processing_a3',\
                        'rx_processing_a4','rx_processing_a5','rx_processing_a6']

    # ancillary
    self.ancArrListL2A=['l2a_alg_count']


    ## L2B
    # array per shot
    self.shotArrListL2B=['algorithmrun_flag','beam','channel','cover','cover_z','delta_time','fhd_normal','l2a_quality_flag','l2b_quality_flag',\
                         'master_frac','master_int','num_detectedmodes','omega','pai','pai_z','pavd_z','pgap_theta','pgap_theta_error','rg',\
                         'rh100','rhog','rhog_error', 'rhov', 'rhov_error','rossg','rv','rx_range_highestreturn','rx_sample_count',\
                         'rx_sample_start_index','selected_l2a_algorithm','selected_mode','selected_mode_flag','selected_rg_algorithm',\
                         'sensitivity','shot_number','stale_return_flag','surface_flag']

    # geolocation. Element per shot
    self.geoArrListL2B=['degrade_flag','delta_time','digital_elevation_model','elev_highestreturn','elev_lowestmode','elevation_bin0',\
                        'elevation_bin0_error','elevation_lastbin','elevation_lastbin_error','height_bin0','height_lastbin',\
                        'lat_highestreturn','lat_lowestmode','latitude_bin0','latitude_bin0_error','latitude_lastbin','latitude_lastbin_error',\
                        'local_beam_azimuth','local_beam_elevation','lon_highestreturn','lon_lowestmode','longitude_bin0','longitude_bin0_error',\
                        'longitude_lastbin','longitude_lastbin_error','shot_number','solar_azimuth','solar_elevation']

    # ancillary
    self.ancArrListL2B=['dz','l2a_alg_count','maxheight_cuttoff','rg_eg_constraint_center_buffer','rg_eg_mpfit_max_func_evals',\
                        'rg_eg_mpfit_maxiters','rg_eg_mpfit_tolerance','signal_search_buff','tx_noise_stddev_multiplier']

    # Fixed number of elements
    self.strArrListL2B=['land_cover_data','rx_processing']


    # L4A files
    # arrays with one element per shot
    self.shotArrListL4A=['agbd', 'agbd_pi_lower', 'agbd_pi_upper', 'agbd_se', 'agbd_t', 'agbd_t_se',\
                         'algorithm_run_flag', 'beam', 'channel', 'degrade_flag', 'delta_time', 'elev_lowestmode', \
                         'l2_quality_flag', 'l4_quality_flag', 'lat_lowestmode', 'lon_lowestmode',\
                         'master_frac', 'master_int', 'predict_stratum', 'predictor_limit_flag', 'response_limit_flag',\
                         'selected_algorithm', 'selected_mode', 'selected_mode_flag', 'sensitivity', 'shot_number',\
                         'solar_elevation', 'surface_flag', 'xvar']

    self.geoArrListL4A=['elev_lowestmode_a1', 'elev_lowestmode_a10', 'elev_lowestmode_a2', 'elev_lowestmode_a3',\
                        'elev_lowestmode_a4', 'elev_lowestmode_a5', 'elev_lowestmode_a6', 'lat_lowestmode_a1',\
                        'lat_lowestmode_a10', 'lat_lowestmode_a2', 'lat_lowestmode_a3', 'lat_lowestmode_a4',\
                        'lat_lowestmode_a5', 'lat_lowestmode_a6', 'lon_lowestmode_a1', 'lon_lowestmode_a10',\
                        'lon_lowestmode_a2', 'lon_lowestmode_a3', 'lon_lowestmode_a4', 'lon_lowestmode_a5',\
                        'lon_lowestmode_a6', 'sensitivity_a1', 'sensitivity_a10', 'sensitivity_a2', 'sensitivity_a3',\
                        'sensitivity_a4', 'sensitivity_a5', 'sensitivity_a6', 'shot_number', 'stale_return_flag']

    # ancillary
    self.ancArrListL4A=['pft_lut', 'region_lut']

    # model data
    # this needs adding 'model_data'

    # Fixed number of elements
    self.strArrListL4A=['land_cover_data']

    # agbd bits
    self.agbdArrListL4A=['agbd_a1', 'agbd_a10', 'agbd_a2', 'agbd_a3', 'agbd_a4', 'agbd_a5',\
                         'agbd_a6', 'agbd_pi_lower_a1', 'agbd_pi_lower_a10',\
                         'agbd_pi_lower_a2', 'agbd_pi_lower_a3', 'agbd_pi_lower_a4',\
                         'agbd_pi_lower_a5', 'agbd_pi_lower_a6', 'agbd_pi_upper_a1',\
                         'agbd_pi_upper_a10', 'agbd_pi_upper_a2', 'agbd_pi_upper_a3',\
                         'agbd_pi_upper_a4', 'agbd_pi_upper_a5', 'agbd_pi_upper_a6',\
                         'agbd_se_a1', 'agbd_se_a10', 'agbd_se_a2', 'agbd_se_a3',\
                         'agbd_se_a4', 'agbd_se_a5', 'agbd_se_a6', 'agbd_t_a1',\
                         'agbd_t_a10', 'agbd_t_a2', 'agbd_t_a3', 'agbd_t_a4', 'agbd_t_a5',\
                         'agbd_t_a6', 'agbd_t_pi_lower_a1', 'agbd_t_pi_lower_a10',\
                         'agbd_t_pi_lower_a2', 'agbd_t_pi_lower_a3', 'agbd_t_pi_lower_a4',\
                         'agbd_t_pi_lower_a5', 'agbd_t_pi_lower_a6', 'agbd_t_pi_upper_a1',\
                         'agbd_t_pi_upper_a10', 'agbd_t_pi_upper_a2', 'agbd_t_pi_upper_a3',\
                         'agbd_t_pi_upper_a4', 'agbd_t_pi_upper_a5', 'agbd_t_pi_upper_a6',\
                         'agbd_t_se_a1', 'agbd_t_se_a10', 'agbd_t_se_a2', 'agbd_t_se_a3',\
                         'agbd_t_se_a4', 'agbd_t_se_a5', 'agbd_t_se_a6',\
                         'algorithm_run_flag_a1', 'algorithm_run_flag_a10',\
                         'algorithm_run_flag_a2', 'algorithm_run_flag_a3',\
                         'algorithm_run_flag_a4', 'algorithm_run_flag_a5',\
                         'algorithm_run_flag_a6', 'l2_quality_flag_a1',\
                         'l2_quality_flag_a10', 'l2_quality_flag_a2', 'l2_quality_flag_a3',\
                         'l2_quality_flag_a4', 'l2_quality_flag_a5', 'l2_quality_flag_a6',\
                         'l4_quality_flag_a1', 'l4_quality_flag_a10', 'l4_quality_flag_a2',\
                         'l4_quality_flag_a3', 'l4_quality_flag_a4', 'l4_quality_flag_a5',\
                         'l4_quality_flag_a6', 'predictor_limit_flag_a1',\
                         'predictor_limit_flag_a10', 'predictor_limit_flag_a2',\
                         'predictor_limit_flag_a3', 'predictor_limit_flag_a4',\
                         'predictor_limit_flag_a5', 'predictor_limit_flag_a6',\
                         'response_limit_flag_a1', 'response_limit_flag_a10',\
                         'response_limit_flag_a2', 'response_limit_flag_a3',\
                         'response_limit_flag_a4', 'response_limit_flag_a5',\
                         'response_limit_flag_a6', 'selected_mode_a1', 'selected_mode_a10',\
                         'selected_mode_a2', 'selected_mode_a3', 'selected_mode_a4',\
                         'selected_mode_a5', 'selected_mode_a6', 'selected_mode_flag_a1',\
                         'selected_mode_flag_a10', 'selected_mode_flag_a2',\
                         'selected_mode_flag_a3', 'selected_mode_flag_a4',\
                         'selected_mode_flag_a5', 'selected_mode_flag_a6', 'shot_number',\
                         'xvar_a1', 'xvar_a10', 'xvar_a2', 'xvar_a3', 'xvar_a4', 'xvar_a5',\
                         'xvar_a6']


    return

  ###########################################

  def findBounds(self,meanN,stdev,i):
    '''Find the signal start and end'''
    thresh=3.5*stdev+meanN
    minWidth=3
    binList=np.where(self.wave[i]>thresh)
    buff=15

    topBin=0
    for j in range(0,len(binList[0])):
      if (binList[0][j]==(binList[0][j-1]+1))&(binList[0][j]==(binList[0][j-2]+2)):
        topBin=binList[0][j]
        break

    botBin=binList[len(binList)-1]
    for j in range(len(binList[0])-1,0,-1):
      if (binList[0][j]==(binList[0][j-1]+1))&(binList[0][j]==(binList[0][j-2]+2)):
        botBin=binList[0][j]
        break

    return(self.z[botBin]-buff,self.z[topBin]+buff)

  ###########################################

  def writeCoords(self):
    for i in range(0,len(self.lon)):
      print(self.lon[i],self.lat[i])


# end of gediData class
###########################################


###########################################
# read the command line

if __name__ == '__main__':
  def gediCommands():
    '''
    Read commandline arguments
    '''
    p = argparse.ArgumentParser(description=("Writes out properties of GEDI waveform files"))
    p.add_argument("--input",dest="inName",type=str,help=("Input GEDI HDF5 filename"))
    p.add_argument("--bounds", dest ="bounds", type=float,nargs=4,default=[-100000000,-100000000,100000000000,10000000000], help=("Bounds to plot between. minX minY maxX maxY"))
    p.add_argument("--output",dest="output",type=str,default='teast.h5',help=("Output filename"))
    cmdargs = p.parse_args()
    return cmdargs


###########################################
# the main block

if __name__ == '__main__':
  # read the command line
  cmdargs=gediCommands()
  inName=cmdargs.inName
  bounds=cmdargs.bounds
  output=cmdargs.output

  # read data
  gedi=gediData(filename=inName,minX=bounds[0],maxX=bounds[2],minY=bounds[1],maxY=bounds[3],outName=output)

