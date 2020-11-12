## Stochastic Gravitational Wave Background Probes combination (SGWBProbe)

This Python 3 package, accompanying the paper Campeti, Komatsu, Poletti & Baccigalupi 2020 (https://arxiv.org/abs/2007.04241), allows to compute the binned sensitivity curves and error bars on the fractional
energy density of gravitational waves in the Universe at the present time, <img src="https://render.githubusercontent.com/render/math?math=\Omega_{GW} h^2">, as a function of frequency for a wide selection of experiments, including CMB B-mode experiments, laser and atomic interferometers and Pulsar Timing Arrays (PTA).  

The list of available experiments currently includes:
- LiteBIRD 
- Laser Interferometer Space Antenna (LISA)
- Einstein Telescope (ET)
- Square Kilometer Array (SKA)
- Big Bang Observer (BBO)
- Deci-hertz Interferometer Gravitational wave Observatory (DECIGO)
- <img src="https://render.githubusercontent.com/render/math?math=\mu">Ares
- Decihertz Observatory (DO),
- Atomic Experiment for Dark Matter and Gravity Exploration in Space (AEDGE)

The forecasts are presented for two representative classes of models for the stochastic background of primordial GW:
- the quantum vacuum fluctuation in the metric from single-field slow-roll inflation 
- the source-induced tensor perturbation from the spectator axion-SU(2) inflation models. 

### Scripts reproducing paper's figures 
We provide python scripts that can be used to quickly reproduce each the figures found in the paper Campeti, Komatsu, Poletti & Baccigalupi 2020 (https://arxiv.org/abs/2007.04241), that is:

- [Fig_1.py](https://github.com/pcampeti/SGWBProbe/blob/main/Fig_1.py)
- [Fig_2.py](https://github.com/pcampeti/SGWBProbe/blob/main/Fig_2.py)
- [Fig_3.py](https://github.com/pcampeti/SGWBProbe/blob/main/Fig_3.py)
- [Fig_4-7.py](https://github.com/pcampeti/SGWBProbe/blob/main/Fig_4-7.py)
- [Fig_8.py](https://github.com/pcampeti/SGWBProbe/blob/main/Fig_8.py)
- [Fig_9-13.py](https://github.com/pcampeti/SGWBProbe/blob/main/Fig_9-13.py)
- [Fig_14-16.py](https://github.com/pcampeti/SGWBProbe/blob/main/Fig_14-16.py)
- [Fig_17.py](https://github.com/pcampeti/SGWBProbe/blob/main/Fig_17.py)
- [Fig_18-20.py](https://github.com/pcampeti/SGWBProbe/blob/main/Fig_18-20.py)
- [Fig_21.py](https://github.com/pcampeti/SGWBProbe/blob/main/Fig_21.py)
- [Fig_23.py](https://github.com/pcampeti/SGWBProbe/blob/main/Fig_23.py)

where the name of each script indicates the respective number of the figures in the paper.
The figures produced by these scripts are collected in the folder `figures/`


### Available sensitivity curves and response functions
The strain sensitivity curves as a function of frequency for the experiments used in the paper are included in the `files/` folder in `.npz` format.
Their content can be easily unpacked through a simple script, e.g for the LISA experiment 
```python
LISA = np.load(op.join(op.dirname(__file__),'files/S_h_LISA_xcosmo.npz'))
LISA_freq = LISA['x']
LISA_strain = LISA['y']

```
and used in your own Python code.

The folder `files/Responses/` contains the response function as a function of frequency in `.npy` format for the laser interferometers computed as described in Appendix A of the paper. Their content can be easily unpacked through a simple script, e.g for the LISA experiment
```python
Resp_LISA = np.load(op.join(op.dirname(__file__),'files/Responses/Resp_LISA.npy'))
freq_LISA = np.load(op.join(op.dirname(__file__),'files/Responses/f_R_LISA.npy'))

```
and used in your own Python code.


The error bars for CMB experiments are computed using a Fisher matrix approach. The Fisher matrices (computed as described in Section 3 of the paper) are gathered in the folder `files/LiteBIRD_Fisher_matrices/`. 
The name of each file contains the binning width and the specific power spectrum model used to compute them (e.g. the file `Fisher_1.3_r0.npy` contains the Fisher matrix for a binning width <img src="https://render.githubusercontent.com/render/math?math=\Delta\ln k = 1.3"> for the single-field slow-roll model with <img src="https://render.githubusercontent.com/render/math?math=\r=0">).


### The `sgwbprobe` package
The `sgwbprobe` package provided here contains the modules
- `SGWB_Signal.py` containing the class `Signal_GW` needed to compute the energy density 
for the SU(2) Axion model of inflation and the standard signle-field slow-roll 
one. More details on the input parameters for the `SGWB_Signal` class are given below.
- `Binned_errors.py`containing the class `Binned_GW` needed to compute the binned 
sensitivity given the sensitivity curve for a GW experiment. More details on the input parameters for the `Binned_GW` class are given below.
- `error_boxes.py` containing the function `make_error_boxes` used to plot the error 
rectangles.
- `foregrounds.py` containing the functions used to compute the analytical 
approximations for the GWD, EGWD, MBHB and BBH+BNS foreground component for the
interferometers and PTA (see Section 4.2.1 in the paper).
- `effective_degrees.py` containing the function `g_of_k` used to compute the number of 
effective degrees of freedom <img src="https://render.githubusercontent.com/render/math?math=g_{*\rho}"> and <img src="https://render.githubusercontent.com/render/math?math=g_{*s}"> as a function of 
wavenumber k.
- `interpolate.py` containing the functions `log_interp1d` and `log_interp1d_fill`, which can 
be used to interpolate in logarithimc space.

#### The `SGWB_Signal` class
The `SGWB_Signal` class contains methods useful to compute the energy density of gravitational 
waves for the single-field slow-roll and axion-SU(2) models described in 
Section 2 of the paper.

 The input parameters are:   
 - `r_vac`: float. Tensor-to-scalar ratio for quantum vacuum fluctuations (named simply r
        in the paper).
        
 - `n_T`: float (optional). Tensor spectral index. If None is calculated from the
        inflationary consistency relation.
        
 - `r_star`: float (optional). Parameter of the axion-SU(2) spectrum (see Sec.2.2).   
        
 - `k_p`: float (optional). Parameter of the axion-SU(2) spectrum (see Sec.2.2).
        
 - `sigma`: float (optional). Parameter of the axion-SU(2) spectrum (see Sec.2.2).
        
 - `axion`: Boolean (optional), defaults to `None`. If `True` computes Omega_GW for the axion-SU(2) model, otherwise for the
        standard single-field slow-roll model.
        
 - `k_star`: float (optional). Pivot scale of the tensor power spectrum. Default value 0.05.
        
 - `running`: Boolean (optional). If `True` includes the running in the tensor power spectrum according to
        the inflatonary consistency relation (see Sec.2.1).

#### The `Binned_GW` class
The `Binned_GW` class contains of methods used to compute the binned sensitivity curve given the sensitivity curve of a GW observatory.

The input parameters are:
- `name_exp`: string. 
            Name of the current experiment. 
            
- `kmin`: float.
            Minimum wavenumber k of the signal  
            
- `k`: numpy array.
            array of wavenumbers for the signal.
            
- `N_bins`: integer.
            Number of bins in which the signal is binned along the whole k 
            range, only for plotting purposes. `N_bins=80` is usually enough to
            cover the range from 10^-18 to 10^4 Hz, you should increase it if
            the you are using smaller bins or decrease it for larger bins.
            
- `delta_log_k`: float, the logarithm bin step.
        
- `omega_gw`: numpy array.
            The GW signal (gravitational wave energy density in GWs today in 
            the Universe) computed from the Signal_GW class.
            
- `kmin_sens`: float.
            The minimum k of the sensitivity range of the experiment.
            
- `N_bins_sens`: integer. 
            The number of bins in the experiment sensitivity range.
            
- `T_obs`: float. 
            The mission observation time in seconds.
            
- `tensor_spect`: numpy array (optional). 
            Input tensor power spectrum for the model analyzed (Eqs.(2.1) 
            and (2.5)).
            
- `k_sens`: numpy array.
            Array of wavenumbers k representing the band-width of a 
            given experiment.
            
- `sens_curve`: numpy  array. 
            The experiment instrumental strain sensitivity S_h from Eq.(4.13). 
            
- `CMB`: Boolean (optional). 
            `True` if we are computing a CMB sensitivity curve, `False` otherwise.
            
- `F`: numpy ndarray with shape (N_bins_sens, N_bins_sens) (optional). 
            The CMB sensitivity Fisher matrix (optional, necessary only if we 
            want to compute the sensitivity for a CMB experiment).
            
- `A_S`: float (optional).
            Amplitude of the scalar perturbations spectrum, defaults to 
            2.1e-9.
            
- `interp`: Boolean (optional). 
            `True` if you want to interpolate along the instrument bandwidth the 
            input instrumental strain sensitivity 
            curve (`sens_curve`).
            
- `n_det` : integer. 
            Number of detectors in the cross-correlation for interferometers 
            (see Eq.(4.10)), defaults to 1.
            
- `fgs`: Boolean (optional). 
            `True` if you want error bars including foregrounds residuals.
            
- `sigma_L`: float. 
            Fractional uncertainty on the amplitude of the BBH+BNS foreground 
            given by an external experiment (see Sec.4.2.2).
            
- `cosmic_var`: Boolean (optional). 
            `True` if you want to include the cosmic variance in the 
            interferometer SNR (see Ref.[5]).
            
- `f_R`: numpy array (optional).
            Frequency band for the response function necessary only if you want
            also the cosmic variance for the interferometers. 
            
- `R_auto`: numpy array.
            The frequency response <img src="https://render.githubusercontent.com/render/math?math=R_{II}"> for the interferometer.
            
- `R_12`: numpy array. The frequency response <img src="https://render.githubusercontent.com/render/math?math=R_{IJ}"> with <img src="https://render.githubusercontent.com/render/math?math=I<J"> for the interferometer.


### Example of usage: binned error bars plots
We first need to load and unpack the instrumental strain sensivity curves as a function of frequency, e.g. for LISA
```python
LISA = np.load(op.join(op.dirname(__file__),'files/S_h_LISA_xcosmo.npz'))
LISA_freq = LISA['x']
LISA_strain = LISA['y']

```
We then choose a value for the data-taking efficiency and mission observation time for the experiment, 
```python
eff_LISA = 0.75
LISA_T_obs = 4 * year_sec * eff_LISA

```
We create an instance of the `Signal_GW` class, e.g. for an axion-SU(2) with the following parameters (called AX1 model in the paper)
```python
class_axion1 = Signal_GW(r_vac=1e-5, r_star=400, k_p=1e13, sigma=8.1, axion=True, running=True)

```
and choose the wavenumber range over which we want to plot the spectrum of GWs, e.g. <img src="https://render.githubusercontent.com/render/math?math=10^{-5} - 10^{20}\,Mpc^{-1}">  
```python
k = np.logspace(np.log10(1e-5), np.log10(1e20), 100000)

```
in order to generate the primordial signal <img src="https://render.githubusercontent.com/render/math?math=\Omega_{GW} h^2"> we are interested in
```python
omega_gw = class_axion1.analytic_omega_WK(k)

```
We can now instantiate the class `Binned_GW` with the chosen specifications and signal for LISA:
```python
sens_curve_LISA = np.array(LISA_strain)
k_LISA = np.array(LISA_freq) * 6.5e14

class_binned = Binned_GW(
                         name_exp = 'LISA',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.2,
                         sens_curve=sens_curve_LISA,
                         omega_gw=omega_gw,
                         k_sens=k_LISA,
                         kmin_sens=1.21303790e+10,
                         N_bins_sens=7,
                         T_obs=LISA_xcosmo_T_obs,
                         interp=True,
                         n_det = 1.,
                         sigma_L=1.0
                         )

```
in this case without foreground contamination.
Now, we can obtain the signal binned in frequency from the method `Omega_GW_binning` 
```python
binned_signal_y, binned_signal_x = class_binned.Omega_GW_binning()

```
which we will then plot, and get also the inputs for the `make_error_boxes` function which we will later use to create and plot the error bars, from the `sens_curve_binning` method:
```python
xerr, yerr, bins_mean_point, binned_signal, binned_curve = class_binned.sens_curve_binning()

```
Similarly, we can instantiate the class `Binned_GW` again for the LISA experiments, but adding this time the residual foreground contamination:
```python
class_binned_fgs = Binned_GW(
                         name_exp='LISA_with_fgs',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.2,
                         sens_curve=sens_curve_LISA,
                         omega_gw=omega_gw,
                         k_sens=k_LISA,
                         kmin_sens=1.21303790e+10,
                         N_bins_sens=7,
                         T_obs=LISA_xcosmo_T_obs,
                         n_det = 1.,
                         fgs=True,
                         interp=True,
                         sigma_L=0.1
                         )

xerr_fgs, yerr_fgs, bins_mean_point_fgs, binned_signal_fgs, binned_curve_fgs = class_binned_fgs.sens_curve_binning()

```
Finally we can plot on the same figure the binned signal,
```python
plt.loglog(np.array(binned_signal_x)/6.5e14, binned_signal_y, color='blue',label=r'Axion Signal $r_{\star}=400$, $k_{p}=10^{15}$ $Mpc^{-1}$, $\sigma=9.1$',linewidth=1.0, zorder=18)

```
the error bars for the foreground-less case
```python
_ = make_error_boxes(ax, np.array(bins_mean_point)/6.5e14, binned_signal, xerr/6.5e14, yerr, facecolor='b', alpha=0.7, zorder=10)
```
and the error bars including foreground residuals
```python
_ = make_error_boxes(ax, np.array(bins_mean_point_fgs)/6.5e14, binned_signal_fgs, xerr_fgs/6.5e14, yerr_fgs, facecolor='b', alpha=0.4, zorder=9)

```
Similarly, we can obtain error bars for CMB experiments: in this case the class `Binned_GW` needs to be instantiated with the appropriate Fisher matrix and the tensor power spectrum obtained from the `total_spect` method of the `Signal_GW` class
```python
Fisher_axion = np.load(op.join(op.dirname(__file__),'files/LiteBIRD_Fisher_matrices/Fisher_1.2_AX1.npy'))  
power_spectrum_axion = class_axion1.total_spect(k)

class_binned_axion_CMB = Binned_GW(
                         name_exp = 'LiteBIRD',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.2,
                         omega_gw=omega_gw,
                         kmin_sens=1e-4,
                         N_bins_sens=6,
                         CMB=True,
                         F=Fisher_axion,
                         tensor_spect=power_spectrum_axion,
                         sigma_L=1.0
                         )                 
                      
xerr_axion, yerr_axion, bins_mean_point_axion, binned_signal_axion, binned_curve_axion = class_binned_axion_CMB.sens_curve_binning()

```

### Example of usage: binned sensitivity curves
A completely similar procedure can be used to generate binned sensitivity curves as in Fig.8 in our paper. In this case we will just use the binned <img src="https://render.githubusercontent.com/render/math?math=\Omega_{GW} h^2"> sensitivity curve also given by the `sens_curve_binning` method:
```python
ax.loglog(np.array(bins_mean_point)/6.5e14, binned_curve[:len(bins_mean_point)], label='LISA w/o fgs', linewidth='1.5', color=sns.xkcd_rgb["black"], linestyle='--')

ax.loglog(np.array(bins_mean_point_fgs)/6.5e14, binned_curve_fgs[:len(bins_mean_point_fgs)], label='LISA w/ fgs', linestyle='-', linewidth='1.5', color=sns.xkcd_rgb["black"])

```

