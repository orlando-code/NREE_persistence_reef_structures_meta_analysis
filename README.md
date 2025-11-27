This repository documents the code necessary to run the meta-analysis presented in '*Persistence of coral reef structures into the 21st Century*' in Nature Reviews Earth & Environment (202X).


# How to use this repository

If you're in a rush, please check out the `meta-analysis` notebook. This exhibits all the code run from start to finish generating the plots and values detailed in the paper (predominantly in the supplementary material).

If you're interested in running the analysis yourself, read on!

# How to run the analysis

## 0. <ins>Create</ins> a virtual environment with the necessary packages


To get started with the analysis, you need to set up a Python environment with the correct dependencies.
We recommend using **either** Python virtual environments (`venv`) **or** `conda` environments. All required dependencies are specified in `pyproject.toml`.

#### **Option 1: Using `venv` (standard Python virtual environment)**
First, make sure you have Python 3.10+ installed:

```bash
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
python -m pip install --upgrade pip
python -m pip install .    # Installs dependencies from pyproject.toml
```

#### **Option 2: Using `conda`**
If you prefer `conda` (such as Anaconda or Miniconda):

```bash
conda create -n calc_meta python=3.10
conda activate calc_meta
pip install .              # Installs dependencies from pyproject.toml
```

> **Note**
> - The command `pip install .` reads your `pyproject.toml` and installs all runtime (and optionally dev) dependencies.
> - If you run into installation issues with packages that need compilation (e.g., `pyarrow`, `scipy`), try installing the relevant dependencies via `conda` first, then run `pip install .` again.
> - For reproducible builds or pinned versions, edit `pyproject.toml` accordingly.

You're now ready to proceed with data download and organisation!


## 1. <ins>Download</ins> the necessary datasets

The paper and its visualisations make use of a number of datasets. Not all are required for every figure. In order of importance (number of figures relying on them) the datasets are:

a. **Raw extracted calcification data**  
   This data was manually extracted from 122 papers following a systematic review. The dataset generated, detailing measured calcification and bioerosion rates alongside the experimental carbonate conditions, is available at xxx<sup>1</sup>.

b. **Hindcast and forecasted environmental conditions**  
   Past and future pH and sea surface temperature conditions from the Coupled Model Intercomparison Project Version 6 (CMIP6) are provided for the sites of interest as `.csv` files in the `data/climatology/` directory in the repository.

c. **MEOW shapefiles**  
   The Marine Ecoregions of the World (MEOW)<sup>2</sup> dataset is used to visualise the distribution of survey sites. This can be downloaded from the [Resource Watch](https://resourcewatch.org/data/explore/Marine-Ecoregions?section=Discover&selectedCollection=&zoom=3&lat=0&lng=0&pitch=0&bearing=0&basemap=dark&labels=light&layers=%255B%257B%2522dataset%2522%253A%252236803484-c413-49a9-abe2-2286ee99b624%2522%252C%2522opacity%2522%253A1%252C%2522layer%2522%253A%25222dd860af-21be-47c6-8e1d-0b8eb63bfa46%2522%257D%255D&aoi=&page=1&sort=most-viewed&sortDirection=-1) 'Marine Ecoregions' page.

d. **UNEP-WCMC GDCR**  
   The Global Distribution of Coral Reefs from UNEP-WCMC<sup>3</sup> is used to visualise the distribution of reef cover in relation to the survey sites. The distribution of warm-water corals can be downloaded from [https://doi.org/10.34892/t2wk-5t34](https://doi.org/10.34892/t2wk-5t34).

e. **Forecasted coral cover**  
   The response of coral cover to changing environmental conditions from DOI xxxx<<sup>4</sup>> is included as an `.xlsx` file. This is used to contextualise the predicted decline in coral calcification rates in Figure 2 of the main manuscript.

f. **Emissions associated with RCPs**  
   To provide an alternative reference level for calcification rates, Figure 2 can be plotted against CO$_2$ emissions provided in [10.5194/gmd-13-3571-2020](10.5194/gmd-13-3571-2020)<sup>5</sup>.

## 2. <ins>Organise</ins> the datasets into the correct structure. 
Please see below for the necessary repository structure:

```
calcification_meta_analysis/
├── analysis/              # Statistical meta-analysis functions
│   ├── analysis.py        # Core meta-analysis calculations (effect sizes, Cook's distance, etc.)
│   ├── analysis_utils.py  # Utility functions for meta-analysis
│   ├── meta_regression.py # Meta-regression model fitting and predictions
│   ├── metafor.py         # Wrapper for R metafor package integration
│   └── rnative/           # Native R scripts for advanced statistical analyses
│
├── plotting/              # Visualization functions for figures
│   ├── analysis.py        # Plots for meta-regression results and model diagnostics
│   ├── climatology.py     # Climate scenario and climatology visualizations
│   ├── data_exploration.py # Exploratory data analysis plots
│   ├── locations.py       # Geographic and site location maps
│   ├── plot_config.py     # Configuration constants (colors, labels, etc.)
│   └── plot_utils.py      # Shared plotting utilities and helpers
│
├── processing/            # Data cleaning, transformation, and preparation
│   ├── cleaning.py        # Data cleaning and quality control functions
│   ├── processing.py      # Main data processing pipeline
│   ├── carbonate_processing.py # Carbonate chemistry calculations
│   ├── climatology.py     # Climate data processing and merging
│   ├── groups_processing.py # Core grouping and taxonomic processing
│   ├── locations.py       # Location data processing and coordinate assignment
│   ├── taxonomy.py        # Taxonomic classification and mapping
│   └── units.py           # Unit conversion and standardization
│
└── utils/                 # Configuration and utility functions
   ├── config.py          # Project configuration and path settings
    ├── file_ops.py        # File I/O utilities
    ├── r_context_handler.py # R environment management and package imports
    └── utils.py           # General utility functions
|
└── data/
    └── climatology/
        ├── MEOW/            # (c.) currently empty, populate as detailed above
        ├── coral_cover_ph_scenarios_output_table_site_locations.csv    # (b.)
        ├── coral_cover_sst_scenarios_output_table_site_locations.csv   # (b.)
        ├── ph_scenarios_output_table_site_locations.csv                # (b.)
        ├── sst_scenarios_output_table_site_locations.csv               # (b.)
        ├── SUPPLEMENT_DataTables_Meinshausen_6May2020.xlsx             # (f.)
        └── directory_description.md
    ├── UNEP_WCMC/          # (d.) currently empty, populate as detailed above
    ├── coral_cover_locs.xlsx           # (e.)
    ├── extracted_bioerosion_data.xlsx  # (a.)
    ├── extracted_bioerosion_data.xlsx  # (a.)
    └── directory_description.md        # information about files in the directory
└── notebooks/
    └── meta-analysis.ipynb   # runs all processing and plotting code to create paper figures
└── resources/          # (see below for more details)
    ├── api_keys.yaml   # currently missing: create and populate with a Google Maps Geoencoder API key 
    ├── gmaps_locations.yaml   # relevant coordinate pairs for sample locations from Google Maps Geoencoder
    ├── locations.csv   # all relevant coordinate pairs for sample locations (csv format)
    ├── locations.yaml  # all relevant coordinate pairs for sample locations (yaml format)
    ├── mapping.yaml    # file containing the cipher for units mapping, spreadsheet columns etc.
    ├── species_mapping.yaml   # species information from WoRMS
    └── directory_description.md   # information about files in the directory
    
```

## 3. <ins>Run</ins> the code via `notebooks/meta-analysis.ipynb`!

This optionally generates the plots and saves them to a local directory.


### Resources directory

While the following files are created automatically by the analysis pipeline, the complete files are included in the repository for ease (and speed) of use.

`api_keys.yaml`

Where sample sites were not detailed via a pair of valid coordinates these were inferred from location information in the text via the [Google Maps Geocoding API](https://developers.google.com/maps/documentation/geocoding). This service requires a Google account with billing information (although the service itself is free) which is connected via an API key. See the '[Get Started](https://developers.google.com/maps/documentation/geocoding/start?_gl=1*1inh8nb*_up*MQ..*_ga*NzcyNzgzMDAzLjE3NjM0MDUxMTA.*_ga_SM8HXJ53K2*czE3NjM0MDUxMTAkbzEkZzAkdDE3NjM0MDUxMTAkajYwJGwwJGgw*_ga_NRWSTWS78N*czE3NjM0MDUxMTAkbzEkZzAkdDE3NjM0MDUxMTAkajYwJGwwJGgw)' page for more information. By default

`gmaps_locations.yaml`

Contains Google Maps geocoded latitude and longitude pairs for sample or study locations, as obtained via the Google Maps API. Used to standardize coordinates for locations missing explicit coordinate data.

`locations.csv`

A CSV file containing all relevant site coordinates and location names, compiling data from various sources. Used to cross-reference or supplement location information with standardized latitude and longitude formats.

`locations.yaml`

A YAML-formatted version of the full set of location information, including coordinates and human-readable location names. Used throughout the codebase for easy merging and referencing of location data.

`mapping.yaml`

Provides mappings for spreadsheet column names, measurement units, and other key-value lookups to standardize data inputs during processing. This file acts as a central cipher for harmonizing the structure of datasets throughout the analysis pipeline.

`species_mapping.yaml`

Contains taxonomic mapping information including genus, species, family, and higher taxa for each organism, based primarily on data retrieved from the World Register of Marine Species (WoRMS) API. Used to ensure consistent taxonomy assignments throughout the analyses.


`species_mapping.yaml`

Taxonomic data (genus, species, family, taxa) was assigned based on the reported species binomial via the World Register of Marine Species (WoRMS) ‘AphiaRecordsByName’ API (https://www.marinespecies.org/rest/)

### Contact

If you spot something odd or otherwise have any feedback on this codebase, please do get in touch via a [GitHub Issue](https://github.com/orlando-code/NREE_persistence_reef_structures_meta_analysis/issues).

## References

<sup>1</sup> DATA REPO

<sup>2</sup> Mark D. Spalding, Helen E. Fox, Gerald R. Allen, Nick Davidson, Zach A. Ferdaña, Max Finlayson, Benjamin S. Halpern, Miguel A. Jorge, Al Lombana, Sara A. Lourie, Kirsten D. Martin, Edmund McManus, Jennifer Molnar, Cheri A. Recchia, James Robertson, Marine Ecoregions of the World: A Bioregionalization of Coastal and Shelf Areas, BioScience, Volume 57, Issue 7, July 2007, Pages 573–583, https://doi.org/10.1641/B570707

<sup>3</sup> UNEP-WCMC, WorldFish Centre, WRI, TNC (2021). Global distribution of coral reefs, compiled from multiple sources including the Millennium Coral Reef Mapping Project. Version 4.1, updated by UNEP-WCMC. Includes contributions from IMaRS-USF and IRD (2005), IMaRS-USF (2005) and Spalding et al. (2001). Cambridge (UK): UN Environment Programme World Conservation Monitoring Centre. Data DOI: https://doi.org/10.34892/t2wk-5t34

<sup>4</sup> Cover paper

<sup>5</sup> Meinshausen, M., Nicholls, Z. R. J., Lewis, J., Gidden, M. J., Vogel, E., Freund, M., Beyerle, U., Gessner, C., Nauels, A., Bauer, N., Canadell, J. G., Daniel, J. S., John, A., Krummel, P. B., Luderer, G., Meinshausen, N., Montzka, S. A., Rayner, P. J., Reimann, S., Smith, S. J., van den Berg, M., Velders, G. J. M., Vollmer, M. K., and Wang, R. H. J.: The shared socio-economic pathway (SSP) greenhouse gas concentrations and their extensions to 2500, Geosci. Model Dev., 13, 3571–3605, https://doi.org/10.5194/gmd-13-3571-2020, 2020.
