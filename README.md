# SDSS Galaxy Data Collector

A Python tool for collecting and processing galaxy data from the Sloan Digital Sky Survey (SDSS) Data Release 16. This tool is designed for astronomers and machine learning researchers working with galaxy data.

## Features

- Query SDSS DR16 for galaxy data with customizable filters
- Download high-quality galaxy images in multiple bands (g, r, i)
- Process images with advanced techniques for machine learning
- Create augmented versions of galaxy images
- Save comprehensive metadata for each galaxy
- Built-in rate limiting to respect SDSS API constraints
- Create color-magnitude diagrams
- Grid visualization of multiple galaxies
- Support for known galaxies catalog

## Prerequisites

1. Python 3.8 or higher

## Installation

1. Clone the repository:
```bash
git clone https://github.com/plugtheo/astro_data_scaper.git
```

2. Create and activate a virtual environment (recommended):
```bash
# On Windows
python -m venv venv
.\venv\Scripts\activate

# On Unix or MacOS
python -m venv venv
source venv/bin/activate
```

3. Install dependencies:
```bash
pip install -r requirements.txt
```

## Usage

```python
from galaxy_data_collector import GalaxyDataCollector

# Initialize the collector
collector = GalaxyDataCollector(output_dir="data")

# Process galaxies from CSV file
# Get the root directory (one level up from src)
root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
csv_file = os.path.join(root_dir, "data", "YOUR_CSV_FILE.csv")
galaxies = collector.process_galaxy_coordinates(csv_file, limit=None)

# Download multiple SDSS images per galaxy with augmentation
metadata = collector.get_sdss_cutout_images_augmented(
    galaxies, n_versions=3, shift_arcsec=1.0, width=256, height=256
)

# Save metadata for images
metadata_path = os.path.join(collector.output_dir, "galaxy_metadata.json")
collector.save_galaxy_metadata(metadata, metadata_path)
```

## Project Structure

```
astro_data_scaper/
├── src/
│   ├── galaxy_data_collector.py
│   └── main.py
├── data/
│   └── images/
├── requirements.txt
├── README.md
└── .env
```

## Data Fields

The following fields are collected for each galaxy:

### Basic Identification
- `objID`: Unique SDSS object identifier (64-bit integer)
- `ra`: Right Ascension in degrees (J2000) - position on the sky
- `dec`: Declination in degrees (J2000) - position on the sky
- `run`, `rerun`, `camcol`, `field`: SDSS imaging run identifiers used for image retrieval
  - `run`: SDSS imaging run number
  - `rerun`: Processing rerun number
  - `camcol`: Camera column (1-6)
  - `field`: Field number within the run
- `obj`: Object number within the field

### Magnitude Measurements
All magnitudes are in the AB system and represent the brightness of the galaxy in different wavelength bands.

#### Model Magnitudes
- `modelMag_u`, `modelMag_g`, `modelMag_r`, `modelMag_i`, `modelMag_z`: Model magnitudes in each band
  - `u`: Ultraviolet (354.3 nm)
  - `g`: Green (477.0 nm)
  - `r`: Red (623.1 nm)
  - `i`: Near-infrared (762.5 nm)
  - `z`: Infrared (913.4 nm)
- `modelMagErr_u`, `modelMagErr_g`, `modelMagErr_r`, `modelMagErr_i`, `modelMagErr_z`: Errors on model magnitudes

#### PSF Magnitudes
- `psfMag_u`, `psfMag_g`, `psfMag_r`, `psfMag_i`, `psfMag_z`: Point Spread Function magnitudes
  - Measures the total flux assuming the object is point-like
  - Useful for stars and compact objects

#### Fiber Magnitudes
- `fiberMag_u`, `fiberMag_g`, `fiberMag_r`, `fiberMag_i`, `fiberMag_z`: Fiber magnitudes
  - Measures flux within a 3-arcsecond diameter aperture
  - Used for spectroscopic observations

#### Composite Model Magnitudes
- `cModelMag_u`, `cModelMag_g`, `cModelMag_r`, `cModelMag_i`, `cModelMag_z`: Composite model magnitudes
  - Best-fit combination of exponential and de Vaucouleurs profiles
  - Most accurate total magnitude measurement

### Morphological Parameters
- `petroRad_r`: Petrosian radius in r-band (arcseconds)
  - Measures galaxy size
  - Radius containing 90% of the total light
- `fracDeV_r`: Fraction of light fit by de Vaucouleurs profile (0-1)
  - 0: Pure exponential disk (spiral galaxy)
  - 1: Pure de Vaucouleurs profile (elliptical galaxy)
- `expRad_r`: Exponential fit radius in r-band (arcseconds)
  - Scale length of the exponential disk component
- `deVRad_r`: de Vaucouleurs fit radius in r-band (arcseconds)
  - Scale length of the bulge component
- `expAB_r`: Exponential fit axis ratio in r-band (0-1)
  - 1: Perfect circle
  - <1: Elliptical shape
- `deVAB_r`: de Vaucouleurs fit axis ratio in r-band (0-1)
- `expPhi_r`: Exponential fit position angle in r-band (degrees)
  - Angle of the major axis from North
- `deVPhi_r`: de Vaucouleurs fit position angle in r-band (degrees)

### Image Position
- `rowc`, `colc`: Row and column coordinates in SDSS image (pixels)
  - Used for precise image cutout positioning

### Quality Flags
- `type`: Object type
  - 3: Galaxy
  - 6: Star
  - Other values indicate different object types
- `class`: Spectral class (if available)
  - G: Galaxy
  - Q: Quasar
  - STAR: Star
- `clean`: Clean photometry flag
  - 1: Clean photometry
  - 0: Potential issues
- `flags`: Bitmask of photometric processing flags
  - Indicates various processing conditions and potential issues
- `flags2`: Additional photometric flags
  - Extended set of processing flags
- `objc_type`: Object classification
  - More detailed classification than type

### Derived Quantities
- `g_r`: g-r color index (modelMag_g - modelMag_r)
  - Redder galaxies have larger values
  - Useful for galaxy classification
- `u_r`: u-r color index (modelMag_u - modelMag_r)
  - Sensitive to star formation and dust

## SDSS Query Details

The galaxy data is collected using the following SDSS DR16 query from SDSS CasJobs:

```sql
SELECT TOP 100000
    g.objID,
    g.ra, g.dec,
    g.run, g.rerun, g.camcol, g.field,
    g.modelMag_u, g.modelMag_g, g.modelMag_r, g.modelMag_i, g.modelMag_z,
    g.modelMagErr_u, g.modelMagErr_g, g.modelMagErr_r, g.modelMagErr_i, g.modelMagErr_z,
    g.psfMag_u, g.psfMag_g, g.psfMag_r, g.psfMag_i, g.psfMag_z,
    g.fiberMag_u, g.fiberMag_g, g.fiberMag_r, g.fiberMag_i, g.fiberMag_z,
    g.cModelMag_u, g.cModelMag_g, g.cModelMag_r, g.cModelMag_i, g.cModelMag_z,
    g.petroRad_r, g.fracDeV_r, g.expRad_r, g.deVRad_r,
    g.expAB_r, g.deVAB_r, g.expPhi_r, g.deVPhi_r,
    g.rowc, g.colc, g.obj,
    g.clean, g.flags
FROM Galaxy AS g
JOIN SpecObj AS s ON s.bestobjid = g.objid
WHERE
    g.petroRad_r > 10  
    AND g.modelMag_u BETWEEN 0 AND 19.6  
    AND g.modelMag_g BETWEEN 0 AND 17    
    AND g.modelMag_r < 18
    AND g.modelMagErr_r < 0.2
```

### Query Explanation

This query selects galaxies from the SDSS DR16 database with the following criteria:

1. **Size Selection**:
   - `g.petroRad_r > 10`: Selects galaxies with a Petrosian radius greater than 10 arcseconds
   - This ensures we get relatively large, well-resolved galaxies

2. **Magnitude Range**:
   - `g.modelMag_u BETWEEN 0 AND 19.6`: Ultraviolet magnitude range
   - `g.modelMag_g BETWEEN 0 AND 17`: Green magnitude range
   - `g.modelMag_r < 18`: Red magnitude cutoff
   - These cuts select galaxies bright enough to be well-measured but not so bright as to be saturated

3. **Quality Control**:
   - `g.modelMagErr_r < 0.2`: Ensures good quality measurements in the r-band
   - The query joins with the `SpecObj` table to ensure we only get galaxies with spectroscopic observations

4. **Data Selection**:
   - The query retrieves all available photometric measurements (model, PSF, fiber, and composite model magnitudes)
   - Includes morphological parameters and quality flags
   - Limited to 100,000 objects to manage data volume

This query is designed to select a high-quality sample of well-resolved galaxies with good photometric measurements and spectroscopic data.

## Application Purpose and Process

This application is designed to collect and process galaxy data from the Sloan Digital Sky Survey (SDSS) for scientific analysis and machine learning applications. Here's how it works:

### 1. Data Collection Process

1. **Initial Query**:
   - The application starts by querying SDSS DR16 for galaxies meeting specific criteria
   - Uses the query described above to select well-resolved, bright galaxies
   - Ensures data quality through magnitude and error cuts
   - Requires spectroscopic observations for additional validation

2. **Image Download**:
   - For each selected galaxy, downloads three-band images (g, r, i)
   - Uses precise coordinates and run/field information for accurate image retrieval
   - Field size is calculated based on galaxy size (petroRad_r)
   - Images are downloaded in FITS format and converted to JPG

3. **Image Processing**:
   - Applies logarithmic scaling to enhance faint features
   - Normalizes each band using percentile-based scaling
   - Creates RGB composite images (i->R, r->G, g->B)
   - Saves high-quality JPG images (95% quality)

### 2. Data Augmentation

The application includes an augmentation feature that:
- Creates multiple versions of each galaxy image
- Applies small random shifts to the center position
- Helps in creating robust machine learning datasets
- Default settings: 3 versions per galaxy, 1 arcsecond max shift

### 3. Output Organization

The processed data is organized into:
- `data/images/`: Contains the processed galaxy images
- `data/sdss_cutouts/`: Contains augmented image versions
- JSON metadata files with all galaxy properties
- Separate metadata for original and augmented images

### 4. Purpose and Applications

This tool is particularly useful for:
- Creating datasets for galaxy classification
- Studying galaxy morphology
- Training machine learning models
- Analyzing galaxy properties
- Creating visualizations of galaxy samples

The combination of high-quality photometric data and processed images makes this tool useful for both scientific research and machine learning applications in astronomy.

## SDSS Data Usage Requirements

When using this project with SDSS data, you must include the following attribution:

### Required Citation
```
Funding for the Sloan Digital Sky Survey IV has been provided by the Alfred P. Sloan Foundation, the U.S. Department of Energy Office of Science, and the Participating Institutions. SDSS-IV acknowledges support and resources from the Center for High-Performance Computing at the University of Utah. The SDSS web site is www.sdss.org.

SDSS-IV is managed by the Astrophysical Research Consortium for the Participating Institutions of the SDSS Collaboration including the Brazilian Participation Group, the Carnegie Institution for Science, Carnegie Mellon University, the Chilean Participation Group, the French Participation Group, Harvard-Smithsonian Center for Astrophysics, Instituto de Astrofísica de Canarias, The Johns Hopkins University, Kavli Institute for the Physics and Mathematics of the Universe (IPMU) / University of Tokyo, the Korean Participation Group, Lawrence Berkeley National Laboratory, Leibniz Institut für Astrophysik Potsdam (AIP), Max-Planck-Institut für Astronomie (MPIA Heidelberg), Max-Planck-Institut für Astrophysik (MPA Garching), Max-Planck-Institut für Extraterrestrische Physik (MPE), National Astronomical Observatories of China, New Mexico State University, New York University, University of Notre Dame, Observatário Nacional / MCTI, The Ohio State University, Pennsylvania State University, Shanghai Astronomical Observatory, United Kingdom Participation Group, Universidad Nacional Autónoma de México, University of Arizona, University of Colorado Boulder, University of Oxford, University of Portsmouth, University of Utah, University of Virginia, University of Washington, University of Wisconsin, Vanderbilt University, and Yale University.
```

For more information about SDSS data usage and licensing, visit [SDSS Data Release License](https://www.sdss.org/collaboration/citing-sdss/).

## License

MIT License 