# Scale factors related to electrons, phothons

## Using the photon scale factor

### Central value
The medium photon ID scale factor for photon pt > 200 GeV is provided in bins of absolute pseudorapidity in `photon_medium_id_sf_v0.root`. The histograms are `photon_medium_id_sf_2017` and `photon_medium_id_sf_2018`. The same central value is used independent of pt.

### Uncertainty
The uncertainty on the SF comes in two parts:

- baseline uncertainty from tag&probe fit variations, etc. This uncertainty is provided as the bin error in the same histograms as the central values above.

- extrapolation unertainty: To take into account residual non-flatness of the SF, an extrapolation uncertainty is applied that increases linearly with photon pt. Like the central values, the slope of the uncertainty is available in the `photon_medium_id_extrap_unc_2017` and `photon_medium_id_extrap_unc_2018` histograms in the same file as before. This uncertainty variation can then we propagated for a given photon by calculating: `SF up / down = SFcentral +/- slope * (photon pt - 150)`.

The two different sources of uncertainty are to be implemented as separate nuisance parameters. They are assumed to be uncorrelated between years.

