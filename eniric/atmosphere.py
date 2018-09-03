import warnings
from os.path import join
from typing import List, Optional

import numpy as np
from astropy import constants as const
from numpy import ndarray

import eniric
import eniric.io_module as io
from eniric.broaden import resolution_convolution
from eniric.utilities import band_limits


class Atmosphere(object):
    """Atmospheric transmission object.

    Stores wavelength and atmospheric transmission arrays.
    Enables telluric masking and accounting for barycentric motion.

    Attributes
    ----------
    wl: ndarray
        Wavelength array
    transmission: ndarray
        Atmospheric transmission (between 0 and 1)
    std: ndarray
        Standard deviation of transmission.
    mask: ndarray
        Transmission mask (1's are kept)
    shifted: bool
        Indicate shifted mask

    Constructors
    ----------
    from_file(atmmodel)
        Read in atmospheric model and prepare.
    from_band(band, bary=False)
        Read in atmospheric model for given band.

    Methods
    -------
    to_file(fname, header, fmt)
        Save the atmospheric model to a txt file.
    at(wave)
        Return the transmission value at the closest points to wave.
    wave_select(wl_min, wl_max)
       Slice Atmosphere between two wavelengths.
    band_select(band)
        Slice Atmosphere to a given band.
    copy()
        Make a copy of atmosphere object.
    mask_transmission(depth)
        Mask the transmission below given depth. e.g. 2%
    bary_shift_mask(rv, consecutive_test)
        Sweep telluric mask symmetrically by rv.
    broaden(resolution, *kwargs)
        Instrument broadening of the atmospheric transmission profile.

    Configuration
    -------------
    Two things can be set for the Atmosphere class in the `config.yaml` file
    The path to atmosphere data
    e.g.
        paths:
            atmmodel: "path/to/atmmodel/directory"
    The name for the atmosphere model .txt file
        atmmodel:
            base: "Average_TAPAS_2014"
    """

    def __init__(self, wavelength, transmission, mask=None, std=None, shifted=False):
        assert len(wavelength) == len(
            transmission
        ), "Wavelength and transmission do not match length."
        self.wl = np.asarray(wavelength)
        self.transmission = np.asarray(transmission)
        if std is None:
            self.std = np.zeros_like(wavelength)
        else:
            self.std = np.asarray(std)
        if mask is None:
            self.mask = np.ones_like(wavelength, dtype=bool)
        else:
            self.mask = np.asarray(mask, dtype=bool)
        self.shifted = shifted

    @classmethod
    def from_file(cls, atmmodel: str):
        """Read in atmospheric model and prepare.

        Alternate constructor for Atmosphere.

        Parameters
        ----------
        atmmodel: str
            Name of atmosphere file.
        """
        wav_atm, flux_atm, std_flux_atm, mask_atm = io.pdread_4col(atmmodel)
        wav_atm = wav_atm / 1e3  # conversion from nanometers to micrometers
        mask_atm = np.array(mask_atm, dtype=bool)
        # We do not use the std from the year atmosphere but need it for compatibility.
        shifted = True if "_bary" in atmmodel else False
        return cls(
            wavelength=wav_atm,
            transmission=flux_atm,
            mask=mask_atm,
            std=std_flux_atm,
            shifted=shifted,
        )

    @classmethod
    def from_band(cls, band: str, bary: bool = False):
        """Read in atmospheric model for given band.

        Alternate constructor for Atmosphere. Base on "base" path in config.yaml.

        Parameters
        ----------
        band: str
            Name of atmosphere file.
        bary: bool
            Barycentric shifted mask.
        """

        extension = "_bary.txt" if bary else ".txt"
        atmmodel = join(
            eniric.paths["atmmodel"],
            "{0}_{1}{2}".format(eniric.atmmodel["base"], band, extension),
        )

        try:
            # Try find the band file
            atm = cls.from_file(atmmodel)
        except IOError:
            warnings.warn(
                """Could not find band file for band {0}.
             It is recommend to create this using
                `split_atmosphere.py -b {0}`
                `bary_shift_atmmodel.py -b {0}`
             Trying to load main atmosphere file for now. (will be slower).""".format(
                    band
                )
            )
            full_model = join(
                eniric.paths["atmmodel"], "{0}.txt".format(eniric.atmmodel["base"])
            )
            atm = cls.from_file(full_model)

            # Shorten to band
            atm = atm.wave_select(*band_limits(band))
            if bary:
                atm.bary_shift_mask(consecutive_test=True)
        return atm

    def to_file(
        self, fname: str, header: Optional[List[str]] = None, fmt: str = "%11.8f"
    ):
        """Save the atmospheric model to a txt file.

        Converts micron back into nanometers to be consistent with from_file().

        Parameters
        ----------
        fname: str
            Name of atmosphere file to save to.
        header:
            Header lines to add.
        fmt: str
             String formatting
        """
        if header is None:
            header = ["# atm_wav(nm)", "atm_flux", "atm_std_flux", "atm_mask"]
        return io.pdwrite_cols(
            fname,
            self.wl * 1000,
            self.transmission,
            self.std,
            self.mask.astype(int),
            header=header,
            float_format=fmt,
        )

    def __getitem__(self, item):
        """Index Atmosphere by returning a Atmosphere with indexed components."""
        return Atmosphere(
            wavelength=self.wl[item],
            transmission=self.transmission[item],
            mask=self.mask[item],
            std=self.std[item],
        )

    def at(self, wave):
        """Return the transmission value at the closest points to wave.

        This assumes that the atmosphere model is
        sampled much higher than the stellar spectra.

        For instance the default has a sampling if 10 compared to 3.
        (instead of interpolation)

        Parameters
        ----------
        wave: ndarray
            Wavelengths at which to return closest atmosphere values.
        """
        # Getting the wav, flux and mask values from the atm model
        # that are the closest to the stellar wav values, see
        # https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
        index_atm = np.searchsorted(self.wl, wave)
        wl_len = len(self.wl)
        # replace indexes outside the array, at the very end, by the value at the very end
        # index_atm = [index if(index < len(wav_atm)) else len(wav_atm)-1 for index in index_atm]
        index_mask = index_atm >= wl_len  # find broken indices
        index_atm[index_mask] = wl_len - 1  # replace with index of end.
        return self[index_atm]

    def wave_select(self, wl_min, wl_max):
        """Slice Atmosphere between two wavelengths."""
        wl_mask = (self.wl < wl_max) & (self.wl > wl_min)
        return self[wl_mask]

    def band_select(self, band):
        """Slice Atmosphere to a given Band."""
        wl_min, wl_max = band_limits(band)
        return self.wave_select(wl_min, wl_max)

    def copy(self):
        """Make a copy of atmosphere object."""
        return Atmosphere(
            wavelength=self.wl.copy(),
            transmission=self.transmission.copy(),
            mask=self.mask.copy(),
            std=self.std.copy(),
        )

    def mask_transmission(self, depth: float = 2.0) -> None:
        """Mask the transmission below given depth. e.g. 2%

        Parameters
        ----------
        depth : float (default = 2.0)
            Telluric line depth percentage to mask out.
            E.g. depth=2 will mask transmission deeper than 2%.

        Updates the mask.
        """
        cutoff = 1 - depth / 100.0
        self.mask = self.transmission >= cutoff

    def bary_shift_mask(self, rv: float = 30.0, consecutive_test: bool = False):
        """Sweep telluric mask symmetrically by rv.

        Parameters
        ----------
        rv: float (default=30 km/s)
            Barycentric RV to extend masks in km/s. (Default=30 km/s)
        consecutive_test: bool (default False)
            Checks for 3 consecutive zeros to mask out transmission.

        """
        if self.shifted:
            warnings.warn(
                "Detected that 'shifted' is already True. "
                "Check that you want to rv extend masks again."
            )
        rv_mps = rv * 1e3  # Convert from km/s into m/s

        shift_amplitudes = self.wl * rv_mps / const.c.value
        # Operate element wise
        blue_shifts = self.wl - shift_amplitudes
        red_shifts = self.wl + shift_amplitudes

        bary_mask = []
        for (blue_wl, red_wl, mask) in zip(blue_shifts, red_shifts, self.mask):
            if mask == 0:
                this_mask_value = False
            else:
                # np.searchsorted is faster then the boolean masking wavelength range
                # It returns index locations to put the min/max doppler-shifted values
                slice_limits = np.searchsorted(self.wl, [blue_wl, red_wl])
                slice_limits = [
                    index if (index < len(self.wl)) else len(self.wl) - 1
                    for index in slice_limits
                ]  # Fix searchsorted end index

                mask_slice = self.mask[slice_limits[0] : slice_limits[1]]

                if consecutive_test:
                    # Mask value False if 3 or more consecutive zeros in slice.
                    len_consec_zeros = consecutive_truths(~mask_slice)
                    if np.all(
                        ~mask_slice
                    ):  # All pixels of slice is zeros (shouldn't get here)
                        this_mask_value = False
                    elif np.max(len_consec_zeros) >= 3:
                        this_mask_value = False
                    else:
                        this_mask_value = True
                        if np.sum(~mask_slice) > 3:
                            print(
                                (
                                    "There were {0}/{1} zeros in this "
                                    "barycentric shift but None were 3 consecutive!"
                                ).format(np.sum(~mask_slice), len(mask_slice))
                            )

                else:
                    this_mask_value = np.bool(
                        np.product(mask_slice)
                    )  # Any 0s will make it 0

                # Checks
                if not this_mask_value:
                    assert np.any(~mask_slice)
                else:
                    if not consecutive_test:
                        assert np.all(mask_slice)
            bary_mask.append(this_mask_value)
        self.mask = np.asarray(bary_mask, dtype=np.bool)
        self.shifted = True

    def broaden(self, resolution: float, fwhm_lim: float = 5, num_procs=None):
        """Instrument broadening of the atmospheric transmission profile.

        This does not change any created masks.

        Parameters
        ----------
        resolution: float
            Instrumental resolution/resolving power
        fwhm_lim: int/float
            Number of FWHM to extend convolution.
        num_procs: Optional[int]
            Number of processors to compute the convolution with. Default = total processors - 1
        """
        self.transmission = resolution_convolution(
            wavelength=self.wl,
            extended_wav=self.wl,
            extended_flux=self.transmission,
            R=resolution,
            fwhm_lim=fwhm_lim,
            num_procs=num_procs,
            normalize=True,
        )


def consecutive_truths(condition: ndarray) -> ndarray:
    """Length of consecutive true values in an bool array.

    Parameters
    ----------
    condition: ndarray of bool
        True False array of a condition.

    Returns
    -------
    len_consecutive: ndarray of ints
        Array of lengths of consecutive true values of the condition.

    Notes
    -----
    Solution found at {http://stackoverflow.com/questions/24342047/
                       count-consecutive-occurences-of-values-varying-in-length-in-a-numpy-array}
    """
    if not np.any(condition):  # No match to condition
        return np.array([0])
    else:
        unequal_consec = np.concatenate(
            ([condition[0]], condition[:-1] != condition[1:], [True])
        )
        where_changes = np.where(unequal_consec)[0]  # indices where condition changes
        len_consecutive = np.diff(where_changes)[
            ::2
        ]  # step through every second to get the "True" lengths.
    return len_consecutive
