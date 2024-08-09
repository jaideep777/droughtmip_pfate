## PlantFATE-SPLASH Contribution to the Drought-MIP Project

This file documents the salient features of the model setup used for the PlantFATE-SPLASH contribution to Drought-MIP.

### Code Version

Simulations are done with code version `dmip_1.1`:

Plant-FATE version:

https://github.com/jaideep777/Plant-FATE/tree/release/drought_mip

SPLASH version: 

https://github.com/jaideep777/rsplash/tree/feature/pf_coupling 

### Key Simulation Features

1. Plant-FATE is coupled with SPLASH.
2. The shortest timescale is daily. Vegetation updates (including GPP, NPP) run every 7 days, while soil water updates run daily.

3. **Representation of diversity**: Diversity is represented in terms of three traits: 
   - Wood density (WD)
   - Max height
   - Xylem P50

   The community-level trait distribution is currently assumed to be unimodal, but with a finite (interspecific) variance. Under high-diversity scenarios, this variance allows the community-level mean (CWM) trait values to change over time following fitness gradients. In principle, this variance also introduces correction terms to the calculation of emergent ecosystem properties, but for the current simulations, we have kept the variance small enough for these corrections to be negligible. The variance only contributes to the evolution of traits by gradually ascending the fitness landscape. 
   Under low-diversity scenarios, traits are held constant at the mode of the aforementioned distribution (variance is set to 0, such that there is no trait evolution). In both scenarios, only the mode appears in the model outputs.

4. Calibration is based primarily on NPP and LAI, and secondarily on AGB, which is not very well predicted. For calibrating against NPP and LAI, the photosynthetic parameter `kphio` (quantum yield efficiency) and the seedling survival parameter `npp_sghalf` (productivity required for 50% seedling survival) are tuned. Except for these two parameters, all other parameters are identical between GYF and TNF.

5. In GYF, best-fit `npp_sghalf` is 0, which means trees are highly shade-tolerant, compared to 0.4 in TNF.

6. Plant-FATE uses two kinds of forcing derived from half-hourly data: 
   - Acclimation forcing, which consists of values calculated a 3-hr means around midday (time of max SW $\pm$ 1.5 hrs), is used to compute acclimated $V_{cmax}$ and $J_{max}$.
   - Instantaneous forcing, which consists of daily means, used to calculate daily photosynthesis and water balance.
   - SPLASH takes the daily means as input. 
   - Plant-FATE computes daytime photosynthesis using daytime mean forcing, computed as 2 * daily means, assuming a 12-hr day, and again downscales the resulting photosynthesis by a factor of 2.


### Key Hydraulic Features:

Hydraulics is better developed in this version compared to the version used for AmazonFACE-MIP. The key features are as follows.

1. Whole-plant / leaf P50, which drives stomatal closure, is coordinated with Xylem P50 as follows:

   \[
   P_{\text{50,leaf}} = \frac{P_{\text{50,xylem}}}{\left(\frac{\log(0.12)}{\log(0.5)}\right)^{\frac{1}{b}}} \approx \frac{P_{\text{50,xylem}}}{3.06}
   \]

   This means there is a cost of foregone photosynthesis if `P50_xylem` is very less negative.

2. Xylem P50 incurs a respiratory cost:

   \[
   \sim c \cdot P_{\text{50,xylem}}^2
   \]

   This means that highly negative `P50_xylem` values are not viable.

3. Whole-plant conductivity is not coordinated with other hydraulic traits owing to lack of empirical support. It is kept constant at \(0.5 \times 10^{-16} \, m^3 m^{-2}\).

4. Mortality depends on diameter, wood density (WD), and loss of xylem conductivity:

   \[
   \mu = \mu_0 + \mu_d + \mu_{\text{hyd}}
   \]

   where:

   - \[
   \mu_0 = m_{\gamma} \left(\frac{\text{WD}}{c_{\text{WD0}}}\right)^{e_{\text{WD}_{\gamma}}}
   \]

   - \[
   \mu_d = c_{D0} \left(\frac{\text{WD}}{c_{\text{WD0}}}\right)^{e_{\text{WD}}} D^{e_{D0}} + c_{D1} \exp\left(-\frac{D}{0.01}\right)
   \]

   - \[
   \mu_{\text{hyd}} = m_{\text{hydraulic}} \cdot \text{atanh}\left(\left(\frac{\psi_{\text{xylem}}}{\psi_{\text{crit xylem}}}\right)^{b_{\text{xylem}}}\right)
   \]

   where `psi_crit_xylem` is the point where 88% of xylem conductivity is lost.

   This mortality function captures the wood-density-mediated growth-mortality tradeoff, such that there is an optimum wood density. Currently, the value of `eWD` is calibrated such that at GYF, the optimum wood density is ~750 g/cc,  

### Model simulation procedure

#### High diversity runs

- First, the HD spinup is performed under HB scenario for ~20000 years until ecosystem properties equlibriate and traits stabilize on fitness maxima. This run stops at 1900 CE. 
- For the HD simulations, the spinup run is continued for 100 years under HB (1900-2000 CE), followed by a CLIMATE run for 10000 years from 2000-12000 CE. 

#### Low diversity runs

- For LD runs, the best adapted trait values are derived as the values averaged over years 1900-2000 CE from the HB-HD run. These traits are kept fixed and a spinup is perfomed for 1000 years until ecosystem properties equilibriate (this ends in 1900 CE).
- The spinup run is continued for 100 years under HB (1900-2000 CE), followed by a CLIMATE run from 2000-2500 CE, with traits held constant throughout. 

- All LD runs use the traits inherited from HD from years 1900-2000. For the LD_HB scenario, we still need the outputs from 1900-2000. These cannot be derived from the HD runs because traits fluctuate on short-timescales in HD, even though they are in long-term equilibrium. So we need to explicitly simulate the LD run between 1900-2000. However, before 1900, the fluxes and size-distribution must have been equilibrated even if traits are fixed. This requires a short (1000 yr) spinup before the LD_HB run. For DY and ND, of course, we continue from HB_LD from year 2000.

### Contents of directories

+ `HB_HD_spin` ==> `HB_20ky_evol_7.0` | 20,000 Year spinup under HB climate with evolution, with trait variances for WD, HMAT, and P50X set to 0.02, 0.02, and 0.1, respectively. 

+ `HB_HD_base` ==> `HB_100y_evol_7.0_cont` | Continuation from HB_HD_spin for 100 years with no change in variances (just for reference, not used)

+ `HB_HD_future` ==> `HB_10ky_evol_7.1_cont_control` | Continuation from `HB_HD_base` for 10,000 years with change in variances (just for control to make sure that change in variances cause unintended changes to mean evolved traits under the same climate, not used)

+ `HB_HD_base_real` ==> `HB_100y_evol_7.1_cont` | Continuation from HB_HD_spin for 100 years with variances set from Amazon-FACE data (WD, HMAT, P50X: 0.03215087, 0.09202125, 0.1672269) 

+ `CL_HD_real` ==> `CL_10ky_evol_7.1_cont` | Continuation from `HB_HD_base_real` for 10000 years with evolution under new climate, with variances set from Amazon-FACE data (WD, HMAT, P50X: 0.03215087, 0.09202125, 0.1672269) 

+ `HB_LD_spin` ==> `HB_1ky_ld_7.1` | Spinup with fixed traits (derived from years 1900-2000 of HD run) for 1000 years.

+ `HB_LD_base` ==> `HB_100y_ld_7.1` | Continuation of `HB_LD_spin` for 100 years

+ `CL_LD` ==> `CL_500y_ld_7.1` | Continuation under new climate from `HB_LD_base` for 500 years

### Observations from model runs

1. Under HB (spinup) runs, LAI and GPP are well predicted, as the model is calibrated to these variables. AGB is somewhat overestimated. 

2. Under HB, ESS wood density is lower in TNF compared to GYF, ESS Max height is higher in TNF, and ESS P50 is more negative in TNF. These trait strategies seem to make sense given that TNF is drier.  

3. Under DY/ND, we see a slight increase in NPP in GYF, but a decrease in TNF. This increase in GYF is likely due to slightly higher average PPFD in these scenarios. LAI does not change appreciably in either site - maybe decreases marginally. ET and volumetric SWC both decrease substantially in both GYF and TNF in the ND scenario, but ET increases slightly in GYF-DY. 

4. Under ND, traits evolve towards lower WD, more negative P50, and somewhat lower max height. The trait change in GYF is general much smaller. 

5. High diversity does lead to a slight recovery of loss in basal area compared to the low-diversity scenario. In GYF-DY, where basal area increases under DY, high-diversity causes a lower increase.

6. However, high diversity does NOT lead to a recovery of AGB loss. On the contrary, the adaptations for higher drought (e.g., lower wood density) lead to an ever lower ESS AGB despite recvoery of basal area. This means that under ND, not only does more aridity lead to a loss of standing biomass, but the community becomes susceptible to invasion by species which further contribute to reduction of biomass / C stock.

7. Plant-FATE substantially underpredicts max height (20 against 35 in GYF, 23 against 40 in TNF), this needs to be investigated.
