import numpy as np
class Sn:

    def __init__(self):

        self.electron_intensity_ratio = 1.137/5.60
        self.electron_amplitude_radio = 390.872/387.461
        self.electron_region = [450.0,850.0]
        self.higher_electron_region = [622.0,700.0]
        self.low_energy_region = [0.0,100.0]
        self.xray_region = [33.0,50.0]

        alpha_xray_amplitudes = [24.002, 24.21]
        alpha_xray_intensities = [28.0, 51.8]
        alpha_xray_amplitude = np.average(alpha_xray_amplitudes,weights = alpha_xray_intensities)
        alpha_xray_intensity = sum(alpha_xray_intensities)
        
        beta_xray_amplitudes = [27.238, 27.276, 27.863]
        beta_xray_intensities = [4.66, 9.0, 2.39]
        beta_xray_amplitude = np.average(beta_xray_amplitudes,weights = beta_xray_intensities)
        beta_xray_intensity = sum(beta_xray_intensities)

        self.xray_intensity_ratio = beta_xray_intensity/alpha_xray_intensity
        self.xray_amplitude_radio = beta_xray_amplitude/alpha_xray_amplitude