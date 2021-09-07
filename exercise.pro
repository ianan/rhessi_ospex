; Find the flare
flares = hsi_select_flare(peak_time_range=['27-Apr-2005 21:00:00','28-Apr-2005 00:00:00'],$
 	peak_countrate_range=[200,500],sflag=1,/structure)

help,flares,/str
;** Structure HSI_FLARELISTDATA_EXT, 26 tags, length=248, data length=237:
;   ID_NUMBER       LONG           5042707
;   START_TIME      DOUBLE      8.3064348e+008
;   END_TIME        DOUBLE      8.3064359e+008
;   PEAK_TIME       DOUBLE      8.3064351e+008
;   BCK_TIME        DOUBLE    Array[2]
;   IMAGE_TIME      DOUBLE    Array[2]
;   ENERGY_RANGE_FOUND
;                   FLOAT     Array[2]
;   ENERGY_HI       FLOAT     Array[2]
;   PEAK_COUNTRATE  FLOAT           285.000
;   BCK_COUNTRATE   FLOAT           11.5867
;   TOTAL_COUNTS    FLOAT           53601.0
;   PEAK_CORRECTION FLOAT           2.22660
;   TOTAL_CORRECTION
;                   FLOAT           2.20253
;   POSITION        FLOAT     Array[2]
;   FILENAME        STRING    '                                                                                '
;   FLAGS           BYTE      Array[32]
;   SEG_INDEX_MASK  BYTE      Array[18]
;   EXTRA_BCK       FLOAT     Array[3, 3]
;   SFLAG1          BYTE         1
;   ACTIVE_REGION   INT              0
;   GOES_CLASS      STRING    '                                                                                '
;   RADIAL_OFFSET   FLOAT          0.000000
;   DURATION        FLOAT           116.000
;   X_POSITION      FLOAT          -677.358
;   Y_POSITION      FLOAT          -37.2008
;   RADIAL_DIST     FLOAT           678.379


; Observing Summary Plot -> Time Range for Spectra
; Start by getting observing summary plot for 10 mins before and after flare times from flare list
timer=[flares.start_time-600,flares.end_time+600]
ltc=hsi_obs_summary(obs_time_interval=timer)
ltc->plotman,/saa,/flare,/night,/corrected,/ylog

; Probably want time range of about
timer='27-Apr-05 '+['22:15','22:30']
ltc->set, obs_time_interval=timer
ltc->plotman,/saa,/flare,/night,/corrected,/ylog

; Make the spectrum and SRM files
os = hsi_spectrum()
os-> set, decimation_correct= 1
os-> set, obs_time_interval= timer
os-> set, pileup_correct= 0
os-> set, seg_index_mask= [1,0,1,1,0,1,0,1,1,0,0,0,0,0,0,0,0,0]
os-> set, sp_data_unit= 'Flux'
os->set, sp_energy_binning=3.+findgen(142)/3.
os-> set, sp_semi_calibrated= 0
os-> set, sp_time_interval= 4
os-> set, sum_flag= 1
os-> set, use_flare_xyoffset= 1
os->filewrite, /fits, /buildsrm,all_simplify=0, /create

; load and do the forward fitting in OSPEX
o=ospex()
; if you want it fully automated
o->set, spex_fit_manual=0, spex_fit_reverse=0, $
	spex_fit_start_method='previous_int'
o->set, spex_autoplot_enable=0, spex_fitcomp_plot_resid=0,$
	spex_fit_progbar=0
o-> set, fit_function='vth+thick2'
o-> set, fit_comp_spectrum= ['full', '']
o-> set, fit_comp_model= ['chianti', '']
o-> set, spex_specfile='hsi_spectrum_20050427_221500.fits'
o-> set, spex_drmfile='hsi_srm_20050427_221500.fits'
o-> set, spex_bk_time_interval=['27-Apr-2005 22:16:00', '27-Apr-2005 22:17:24']
o-> set, spex_bk_order=0
o-> set, spex_fit_time_interval=['27-Apr-2005 22:18:24','27-Apr-2005 22:18:36.000']
o-> set, mcurvefit_itmax= 50
o->set, spex_erange=[3,10]
o->set, fit_comp_free=[1,1,0,0,0,0,0,0,0]
o->set, fit_comp_param=[1e-3,1.5,1,1e-2,6,500,20,15,500]
o->dofit
o->set, spex_erange=[10,40]
o->set, fit_comp_free=[0,0,0,1,1,0,0,1,0]
o->dofit
o->set, spex_erange=[7,40]
o->set, fit_comp_free=[1,1,0,1,1,0,0,1,0]
o->dofit
o->set, spex_erange=[3,40]
o->set, fit_comp_free=[1,1,0,1,1,0,0,0,0]
o->dofit

params=o->get(/spex_summ_params)
print,params
;0.00432079      1.42590      1.00000      2.48007      5.74665      500.000      20.0000      12.8798      500.000

; Thermal Energy
vol=6d26
em=params[0]*1d49
tmk=params[1]*1d6/0.086164
kb=1.38065d-23
j2erg=1d7
energy_th=3*sqrt(em*vol)*kb*tmk*j2erg
; Non-thermal Energy
nume=params[3]*1d35
del=params[4]
ec=params[7]
power_nn=1.6d-9*(del-1)*nume*ec/(del-2)
time_diff=(o->get(/spex_fit_time_interval))[1]-(o->get(/spex_fit_time_interval))[0]
energy_nn=power_nn*time_diff

; Plot the result in count space
o->plot_spectrum,/show_fit,/bksub,spex_units='flux'
; Plot the result in photon spce
o->plot_spectrum,/show_fit,/bksub,/photon,spex_units='flux'



; what about continuum only thermal and gaussian line?
obj-> set, fit_function= 'vth+thick2+line'
obj-> set, fit_comp_spectrum= ['continuum', '', '']
obj-> set, fit_comp_params= [0.00390593, 1.73820, 1, 0.490420, 5.47943, 500, 20, 17.5099, 500, 324.101, 6.87546, 0.1]
obj-> set, fit_comp_free_mask= [1B, 1B, 0B, 1B, 1B, 0B, 0B, 1B, 0B, 1B, 1B, 0B]


; what about only detector 4?
os = hsi_spectrum()
os-> set, decimation_correct= 1
os-> set, obs_time_interval= timer
os-> set, pileup_correct= 0
os-> set, seg_index_mask= [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
os-> set, sp_data_unit= 'Flux'
os-> set, sp_energy_binning= 22
os-> set, sp_semi_calibrated= 0
os-> set, sp_time_interval= 4
os-> set, sum_flag= 1
os-> set, use_flare_xyoffset= 1
os->filewrite, /fits, /buildsrm,all_simplify=0, /create,$
	srmfile='hsi_srm_20050427_221500_D4.fits', specfile='hsi_spectrum_20050427_221500_D4.fits'


; fit with vth+thick2 just detector 4
o=ospex()
o->set, spex_fit_manual=0, spex_fit_reverse=0, $
	spex_fit_start_method='previous_int'
o->set, spex_autoplot_enable=0, spex_fitcomp_plot_resid=0,$
	spex_fit_progbar=0
o-> set, fit_function='vth+thick2'
o-> set, fit_comp_spectrum= ['full', '']
o-> set, fit_comp_model= ['chianti', '']
o-> set, spex_specfile='hsi_spectrum_20050427_221500_D4.fits'
o-> set, spex_drmfile='hsi_srm_20050427_221500_D4.fits'
o-> set, spex_bk_time_interval=['27-Apr-2005 22:16:00', '27-Apr-2005 22:17:24']
o-> set, spex_bk_order=0
o-> set, spex_fit_time_interval=['27-Apr-2005 22:18:24','27-Apr-2005 22:18:36.000']
o-> set, mcurvefit_itmax= 50
o->set, spex_erange=[3,10]
o->set, fit_comp_free=[1,1,0,0,0,0,0,0,0]
o->set, fit_comp_param=[1e-3,1.5,1,1e-2,6,500,20,15,500]
o->dofit
o->set, spex_erange=[10,40]
o->set, fit_comp_free=[0,0,0,1,1,0,0,1,0]
o->dofit
o->set, spex_erange=[7,40]
o->set, fit_comp_free=[1,1,0,1,1,0,0,1,0]
o->dofit
o->set, spex_erange=[3,40]
o->set, fit_comp_free=[1,1,0,1,1,0,0,0,0]
o->dofit

print,o->get(/spex_summ_params)
;0.00605546      1.30572      1.00000      3.24113      5.73656      500.000      20.0000      12.0597      500.000
o->plot_spectrum,/show_fit,/bksub,spex_units='flux'
o->plot_spectrum,/show_fit,/bksub,/photon,spex_units='flux'

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
' don't use below for the example again but instead copy and paste for the above fit????
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; fit with vth+thick2+gauss just detector 4
o=ospex()
o->set, spex_fit_manual=0, spex_fit_reverse=0, $
	spex_fit_start_method='previous_int'
o->set, spex_autoplot_enable=0, spex_fitcomp_plot_resid=0,$
	spex_fit_progbar=0
o-> set, fit_function= 'vth+thick2+line'
o-> set, fit_comp_spectrum= ['continuum', '', '']
o-> set, fit_comp_model= ['chianti', '']
o-> set, spex_specfile='hsi_spectrum_20050427_221500_D4.fits'
o-> set, spex_drmfile='hsi_srm_20050427_221500_D4.fits'
o-> set, spex_bk_time_interval=['27-Apr-2005 22:16:00', '27-Apr-2005 22:17:24']
o-> set, spex_bk_order=0
o-> set, spex_fit_time_interval=['27-Apr-2005 22:18:24','27-Apr-2005 22:18:36.000']
o-> set, mcurvefit_itmax= 50
o->set, spex_erange=[3,10]
o->set, fit_comp_free=[1,1,0,0,0,0,0,0,0,1,1,0]
o->set, fit_comp_param=[1e-3,1.5,1,1e-2,6,500,20,15,500,100,6.7,0.1]
o->dofit
o->set, spex_erange=[10,40]
o->set, fit_comp_free=[0,0,0,1,1,0,0,1,0,0,0,0]
o->dofit
o->set, spex_erange=[7,40]
o->set, fit_comp_free=[1,1,0,1,1,0,0,1,0,1,1,0]
o->dofit
o->set, spex_erange=[3,40]
o->set, fit_comp_free=[1,1,0,1,1,0,0,0,0,1,1,0]
o->dofit

print,o->get(/spex_summ_params)

o->plot_spectrum,/show_fit,/bksub,spex_units='flux'
o->plot_spectrum,/show_fit,/bksub,/photon,spex_units='flux'


;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; load and do the forward fitting in OSPEX
o=ospex()
; if you want it fully automated
o->set, spex_fit_manual=0, spex_fit_reverse=0, $
	spex_fit_start_method='previous_int'
o->set, spex_autoplot_enable=0, spex_fitcomp_plot_resid=0,$
	spex_fit_progbar=0
o-> set, fit_function='vth'
o-> set, fit_comp_spectrum= 'full'
o-> set, fit_comp_model= 'chianti'
o-> set, spex_specfile='hsi_spectrum_20050427_221500.fits'
o-> set, spex_drmfile='hsi_srm_20050427_221500.fits'
o-> set, spex_bk_time_interval=['27-Apr-2005 22:16:00', '27-Apr-2005 22:17:24']
o-> set, spex_bk_order=0
o-> set, spex_fit_time_interval=['27-Apr-2005 22:18:54','27-Apr-2005 22:20:06.000']
o-> set, mcurvefit_itmax= 50
o->set, spex_erange=[3,10]
o->set, fit_comp_free=[1,1,0]
o->set, fit_comp_param=[1e-3,1.5,1]
o->dofit

params=o->get(/spex_summ_params)
print,params
;0.00432079      1.42590      1.00000      2.48007      5.74665      500.000      20.0000      12.8798      500.000

; Thermal Energy
vol=6d26
em=params[0]*1d49
tmk=params[1]*1d6/0.086164
kb=1.38065d-23
j2erg=1d7
energy_th=3*sqrt(em*vol)*kb*tmk*j2erg
; Non-thermal Energy
nume=params[3]*1d35
del=params[4]
ec=params[7]
power_nn=1.6d-9*(del-1)*nume*ec/(del-2)
time_diff=(o->get(/spex_fit_time_interval))[1]-(o->get(/spex_fit_time_interval))[0]
energy_nn=power_nn*time_diff

; Plot the result in count space
o->plot_spectrum,/show_fit,/bksub,spex_units='flux'
; Plot the result in photon spce
o->plot_spectrum,/show_fit,/bksub,/photon,spex_units='flux'




