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
effective degrees of freedom <img src="https://render.githubusercontent.com/render/math?math=g_{*rho}"> and <img src="https://render.githubusercontent.com/render/math?math=g_{*s}"> as a function of 
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
        
 - `axion`: Boolean (optional), defaults to None. If True computes Omega_GW for the axion-SU(2) model, otherwise for the
        standard single-field slow-roll model.
        
 - `k_star`: float (optional). Pivot scale of the tensor power spectrum. Default value 0.05.
        
 - `running`: Boolean type (optional). If True includes the running in the tensor power spectrum according to
        the inflatonary consistency relation (see Sec.2.1).



### The `files/` folder
The strain sensitivity curves as a function of frequency for the experiments used in the paper are included in the folder `files/` in `.npz` format.
Their content can be easily unpacked through a simple script 
```python
ET = np.load(op.join(op.dirname(__file__),'files/S_h_ET.npz'))
ET_freq = ET['x']
ET_strain = ET['y']
```
and used in your Python code.

The error bars for CMB experiments are computed using a Fisher matrix approach. The Fisher matrices (computed as described in Section 3 of the paper) are gathered in the folder `files/LiteBIRD_Fisher_matrices/`. 
The name of each file contains the binning width and the specific power spectrum model used to compute them (e.g. the file `Fisher_1.3_r0.npy` contains the Fisher matrix for a binning width <img src="https://render.githubusercontent.com/render/math?math=\Delta\ln k = 1.3"> for the single-field slow-roll model with <img src="https://render.githubusercontent.com/render/math?math=\r=0">).


