## Stochastic Gravitational Wave Background Probes combination (SGWBProbe)

This package, accompanying the paper Campeti, Komatsu, Poletti & Baccigalupi 2020 (https://arxiv.org/abs/2007.04241), allows to compute the binned sensitivity curves and error bars on the fractional
energy density of gravitational waves in the Universe at the present time, <img src="https://render.githubusercontent.com/render/math?math=\Omega_{GW} h^2">, as a function of frequency for a wide selection of experiments, including at the moment:
- LISA
- 
two representative classes
of models for the stochastic background of primordial GW: the quantum vacuum fluctuation
in the metric from single-field slow-roll inflation, and the source-induced tensor perturbation
from the spectator axion-SU(2) inflation models. 

### Scripts to reproduce paper figures 
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

where the name of each script indicates the respective number of the figures in the paper which can reproduce.
The figures produced by these scripts are collected in the folder `figures/`
### `sgwbprobe` package
The `sgwbprobe` package provided contains:
- j

### `files/` folder
The strain sensitivity curves as a function of frequency for the experiments used in the paper are included in the folder `files/` in `.npz` format.
Their content can be easily unpacked through a simple script 
```python
ET = np.load(op.join(op.dirname(__file__),'files/S_h_ET.npz'))
ET_freq = ET['x']
ET_strain = ET['y']
```
and used in your code.

The error bars for CMB experiments are computed using a Fisher matrix approach. The Fisher matrices (computed as described in Section 3 of the paper) are gathered in the folder `files/LiteBIRD_Fisher_matrices/`. 
The name of each file contains the binning width and the specific power spectrum model used to compute them (e.g. the file `Fisher_1.3_r0.npy` contains the Fisher matrix for a binning width <img src="https://render.githubusercontent.com/render/math?math=\Delta\ln k = 1.3"> for the single-field slow-roll model with <img src="https://render.githubusercontent.com/render/math?math=\r=0">).


