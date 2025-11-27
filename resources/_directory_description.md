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

Contains taxonomic mapping information including genus, species, family, and higher taxa for each organism, based primarily on data retrieved from the World Register of Marine Species (WoRMS) ‘AphiaRecordsByName’ API (https://www.marinespecies.org/rest/). Used to ensure consistent taxonomy assignments throughout the analyses.
