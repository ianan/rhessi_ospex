; write the spectrum and srm files
os = hsi_spectrum()
os-> set, decimation_correct= 1
os-> set, obs_time_interval= '6-Jan-04 '+['06:0','06:35']
os-> set, pileup_correct= 0
os-> set, seg_index_mask= [1,0,1,1,0,1,0,1,1,0,0,0,0,0,0,0,0,0]
os-> set, sp_data_unit= 'Flux'
os-> set, sp_energy_binning= 22
os-> set, sp_semi_calibrated= 0
os-> set, sp_time_interval= 4
os-> set, sum_flag= 1
os-> set, use_flare_xyoffset= 1
os->filewrite, /fits, /buildsrm,all_simplify=0


; Fit of nice A1 time
o=ospex()
; if you want it fully automated
o->set, spex_fit_manual=0, spex_fit_reverse=0, $
	spex_fit_start_method='previous_int'
o->set, spex_autoplot_enable=0, spex_fitcomp_plot_resid=0,$
	spex_fit_progbar=0
o-> set, fit_function='vth+thick2'
o-> set, fit_comp_spectrum= ['full', '']
o-> set, fit_comp_model= ['chianti', '']
o-> set, spex_specfile='hsi_spectrum_20040106_060000.fits'
o-> set, spex_drmfile='hsi_srm_20040106_060000.fits'
o-> set, spex_bk_time_interval=['6-Jan-2004 06:32:00', '6-Jan-2004 06:33:20']
o-> set, spex_bk_order=0
o-> set, spex_fit_time_interval=['6-Jan-2004 06:22:20','6-Jan-2004 06:23:00']
o-> set, mcurvefit_itmax= 50
o->set, spex_erange=[6,20]
o->set, fit_comp_free=[1,1,0,0,0,0,0,0,0]
o->set, fit_comp_param=[1e-3,1.5,1,1e-2,3.5,1000,20,40,1000]
o->dofit
o->set, spex_erange=[20,300]
o->set, fit_comp_free=[0,0,0,1,1,0,0,1,0]
o->dofit
o->set, spex_erange=[6,300]
o->set, fit_comp_free=[1,1,0,1,1,0,0,1,0]
o->dofit

params=o->get(/spex_summ_params)
print,params
;o->plot_spectrum,/show_fit,/bksub,spex_units='flux'
;o->plot_spectrum,/show_fit,/bksub,/photon,spex_units='flux'

vol=2d26
em=params[0]*1d49
tmk=params[1]*1d6/0.086164
kb=1.38065d-23
j2erg=1d7
energy_th=3*sqrt(em*vol)*kb*tmk*j2erg
nume=params[3]*1d35
del=params[4]
ec=params[7]
power_nn=1.6d-9*(del-1)*nume*ec/(del-2)
time_diff=(o->get(/spex_fit_time_interval))[1]-(o->get(/spex_fit_time_interval))[0]
energy_nn=power_nn*time_diff

stop

; second easier event?
timer='26-Jun-02 '+['18:15','19:30']
ltc=hsi_obs_summary(obs_time_interval=timer)
ltc->plotman,/saa,/flare,/night,/corrected,/ylog

os = hsi_spectrum()
os-> set, decimation_correct= 1
os-> set, obs_time_interval= '26-Jun-02 '+['18:40','19:08']
os-> set, pileup_correct= 0
os-> set, seg_index_mask= [1,0,1,1,0,1,0,1,1,0,0,0,0,0,0,0,0,0]
os-> set, sp_data_unit= 'Flux'
; 1/3 keV from 3 to 50 keV
os->set, sp_energy_binning=3.+findgen(142)/3.
os-> set, sp_semi_calibrated= 0
os-> set, sp_time_interval= 4
os-> set, sum_flag= 1
os-> set, use_flare_xyoffset= 1
os->filewrite, /fits, /buildsrm,all_simplify=0
;*****************************************************
o=ospex()
o-> set, fit_function='vth+thick2'
o-> set, fit_comp_spectrum= ['full', '']
o-> set, fit_comp_model= ['chianti', '']
o-> set, spex_specfile='hsi_spectrum_20020626_184000.fits'
o-> set, spex_drmfile='hsi_srm_20020626_184000.fits'
o-> set, spex_bk_time_interval='26-Jun-02 '+['18:43:20','18:45:40']
o-> set, spex_bk_order=0
;o-> set, spex_fit_time_interval='26-Jun-02 '+['18:46:36','18:56:36']
o-> set, spex_fit_time_interval= [['26-Jun-2002 18:46:36.000', $
 '26-Jun-2002 18:47:36.000'], ['26-Jun-2002 18:47:36.000', '26-Jun-2002 18:48:36.000'], $
 ['26-Jun-2002 18:48:36.000', '26-Jun-2002 18:49:36.000'], ['26-Jun-2002 18:49:36.000', $
 '26-Jun-2002 18:50:36.000'], ['26-Jun-2002 18:50:36.000', '26-Jun-2002 18:51:36.000'], $
 ['26-Jun-2002 18:51:36.000', '26-Jun-2002 18:52:36.000'], ['26-Jun-2002 18:52:36.000', $
 '26-Jun-2002 18:53:36.000'], ['26-Jun-2002 18:53:36.000', '26-Jun-2002 18:54:36.000'], $
 ['26-Jun-2002 18:54:36.000', '26-Jun-2002 18:55:36.000'], ['26-Jun-2002 18:55:36.000', $
 '26-Jun-2002 18:56:36.000']]
o->set,spex_fit_start_method='previous_iter'
o->set,spex_fit_init_params_only=1
o->set,spex_fit_manual=0

o->set, spex_erange=[3,8]
o->set, fit_comp_free=[1,1,0,0,0,0,0,0,0]
o->set, fit_comp_param=[5e-4,1.2,1,1e-4,7,200,20,20,200]
o->dofit, spex_intervals_tofit=1
o->set, spex_erange=[9,20]
o->set, fit_comp_free=[0,0,0,1,1,0,0,1,0]
o->dofit, spex_intervals_tofit=1
o->set, spex_erange=[3,20]
o->set, fit_comp_free=[1,1,0,1,1,0,0,0,0]
o->dofit, spex_intervals_tofit=1

o->plot_spectrum,/show_fit,/bksub,spex_units='flux',/overlay_back,/show_err
o->plot_spectrum,/show_fit,/bksub,/photon,spex_units='flux',/overlay_back,/show_err

params=o->get(/spex_summ_params)
print,params

vol=2d26
em=params[0,1]*1d49
tmk=params[1,1]*1d6/0.086164
kb=1.38065d-23
j2erg=1d7
energy_th=3*sqrt(em*vol)*kb*tmk*j2erg
nume=params[3,1]*1d35
del=params[4,1]
ec=params[7,1]
power_nn=1.6d-9*(del-1)*nume*ec/(del-2)
time_diff=(o->get(/spex_fit_time_interval))[1]-(o->get(/spex_fit_time_interval))[0]
energy_nn=power_nn*time_diff

;**********************************************

o=ospex()
o-> set, fit_function='vth+thick2'
o-> set, fit_comp_spectrum= ['full', '']
o-> set, fit_comp_model= ['chianti', '']
o-> set, spex_specfile='hsi_spectrum_20020626_184000.fits'
o-> set, spex_drmfile='hsi_srm_20020626_184000.fits'
o-> set, spex_bk_time_interval='26-Jun-02 '+['18:43:20','18:45:40']
o-> set, spex_bk_order=0
;o-> set, spex_fit_time_interval='26-Jun-02 '+['18:46:36','18:56:36']
o-> set, spex_fit_time_interval= [['26-Jun-2002 18:46:36.000', $
 '26-Jun-2002 18:47:36.000'], ['26-Jun-2002 18:47:36.000', '26-Jun-2002 18:48:36.000'], $
 ['26-Jun-2002 18:48:36.000', '26-Jun-2002 18:49:36.000'], ['26-Jun-2002 18:49:36.000', $
 '26-Jun-2002 18:50:36.000'], ['26-Jun-2002 18:50:36.000', '26-Jun-2002 18:51:36.000'], $
 ['26-Jun-2002 18:51:36.000', '26-Jun-2002 18:52:36.000'], ['26-Jun-2002 18:52:36.000', $
 '26-Jun-2002 18:53:36.000'], ['26-Jun-2002 18:53:36.000', '26-Jun-2002 18:54:36.000'], $
 ['26-Jun-2002 18:54:36.000', '26-Jun-2002 18:55:36.000'], ['26-Jun-2002 18:55:36.000', $
 '26-Jun-2002 18:56:36.000']]
o->set,spex_fit_start_method='previous_iter'
o->set,spex_fit_init_params_only=1
o->set,spex_fit_manual=0

o-> set, mcurvefit_itmax= 50
o->set, spex_erange=[3,8]
o->set, fit_comp_free=[1,1,0,0,0,0,0,0,0]
o->set, fit_comp_param=[5e-4,1.2,1,1e-4,7,200,20,20,200]
o->dofit,/all
o->set, spex_erange=[9,20]
o->set, fit_comp_free=[0,0,0,1,1,0,0,1,0]
o->dofit,/all
o->set, spex_erange=[3,20]
o->set, fit_comp_free=[1,1,0,1,1,0,0,0,0]
o->dofit,/all










