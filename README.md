# LPP-SAR-flood-inventory
Multi-event Sentinel-1 SAR flood inventory for the Lam Phra Phloeng watershed (GEE scripts)
# Multi-event SAR Flood Inventory — Lam Phra Phloeng Watershed

Google Earth Engine (GEE) scripts for the paper:
*Tanachaichoksirikun et al., "Multi-event SAR Flood Inventory in a Reservoir-
regulated Monsoon Basin: Adaptive Otsu Thresholding and InSAR Coherence for the
Lam Phra Phloeng Watershed, Thailand"* (submitted to NHESS).

## Scripts
- `1_LPP_Nangka2020.js` / `2_LPP_Dianmu2021.js` / `3_LPP_Noru2022.js`
  — per-event flood detection (Otsu, total-water) + validation + ablation (C0/C0.5/C1/C2)
- `4_LPP_frequency_LULC.js` — 3-event flood-frequency map + land-cover cross-tabulation
- `5_LPP_A1_coherence.js` — InSAR coherence loss and urban-disturbance analysis

## How to run
1. Open each script in the GEE Code Editor (https://code.earthengine.google.com).
2. Ensure these assets are available in your GEE project:
   - Watershed boundary: `projects/lamphra-pleong-project/assets/Watershed_v2`
   - InSAR coherence (Path 164, HyP3/GAMMA): `nangka_corr_ref2608`, `nangka_corr_co0814`,
     `corr_ref_0921`, `corr_ref_2127`, `noru_corr_ref1628`, `noru_corr_co2810`
3. Run; results (Otsu thresholds, areas, accuracy, confusion counts) print to the Console.

## Data
Sentinel-1/-2, JRC GSW, MERIT Hydro, Copernicus GLO-30 DEM and ESA WorldCover are
open via GEE. InSAR coherence products were generated with ASF HyP3 (GAMMA).
