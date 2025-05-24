"""
This module provides functionality to collect and process galaxy data from SDSS DR16.
"""
import json
import logging
import os
import time
from typing import Dict, List

import numpy as np  # For numerical operations and array handling
import pandas as pd  # For data manipulation and CSV handling
import requests  # For HTTP requests to SDSS cutout service
from astropy import \
    coordinates as coords  # For astronomical coordinate handling
from astropy import units as u  # For handling astronomical units
from astropy.io import fits  # For handling FITS image files
from astropy.table import (Table,  # For handling astronomical data tables
                           vstack)
from astroquery.sdss import SDSS  # Main SDSS API interface
from PIL import Image  # For image processing and saving

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class GalaxyDataCollector:
    """
    A class to collect galaxy images and data from SDSS using astroquery.
    
    This class handles:
    - Querying SDSS database for galaxy data
    - Downloading and processing galaxy images
    - Managing rate limits and API requests
    - Saving processed data and metadata
    """
    
    def __init__(self, output_dir: str = "data"):
        """
        Initialize the GalaxyDataCollector.
        
        Args:
            output_dir (str): Directory to save collected data
        """
        self.output_dir = output_dir
        self._create_output_directory()
        self.last_request_time = 0
        self.min_request_interval = 6  # Minimum seconds between requests to respect SDSS rate limits
        
    def _create_output_directory(self):
        """Create output directory structure if it doesn't exist."""
        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "images"), exist_ok=True)
    
    def _rate_limit(self) -> None:
        """
        Implement rate limiting for SDSS API calls.
        
        This method:
        1. Tracks the time between API calls
        2. Enforces a minimum delay between calls
        3. Logs rate limiting information
        
        Note:
            SDSS API has rate limits to prevent overloading their servers.
            This method ensures we stay within those limits by:
            - Waiting at least 1 second between calls
            - Logging when rate limiting is active
        """
        current_time = time.time()
        if hasattr(self, '_last_api_call'):
            time_since_last_call = current_time - self._last_api_call
            if time_since_last_call < 1.0:  # Minimum 1 second between calls
                sleep_time = 1.0 - time_since_last_call
                logger.debug(f"Rate limiting: sleeping for {sleep_time:.2f} seconds")
                time.sleep(sleep_time)
        self._last_api_call = time.time()

    def query_galaxies(self, limit: int = 100) -> List[Dict]:
        """
        Query galaxies using SDSS API.
        
        This method:
        1. Queries SDSS DR16 database for galaxies in a specific region
        2. Applies quality filters to ensure good data
        3. Processes the results into a standardized format
        
        Args:
            limit (int): Maximum number of galaxies to return
            
        Returns:
            List[Dict]: List of galaxies with their properties
            
        Note:
            The query is centered on NGC 4565 (a well-known edge-on spiral galaxy)
            and searches within a 2.5 arcminute radius.
        """
        try:
            self._rate_limit()  # Ensure we respect SDSS API rate limits
            
            # Query for galaxies in a region of the sky
            # Using coordinates for NGC 4565 - a well-known edge-on spiral galaxy
            center_coords = coords.SkyCoord(ra=189.0866*u.degree, dec=25.9877*u.degree)
            
            logger.info("Executing galaxy query using SDSS")
            # SDSS API query with specific fields
            # Each field is documented with its purpose
            results = SDSS.query_region(
                center_coords,
                radius=2.5*u.arcmin,  # Just under the 3 arcmin limit
                spectro=True,  # Only include objects with spectroscopic data
                data_release=16,  # Use SDSS DR16
                fields=[
                    'objID',  # Unique SDSS object identifier
                    'ra', 'dec',  # Right Ascension and Declination (J2000)
                    'run', 'rerun', 'camcol', 'field',  # Imaging run identifiers (for image retrieval)
                    'modelMag_u', 'modelMag_g', 'modelMag_r', 'modelMag_i', 'modelMag_z',  # Model magnitudes in each band
                    'modelMagErr_u', 'modelMagErr_g', 'modelMagErr_r', 'modelMagErr_i', 'modelMagErr_z',  # Errors on model magnitudes
                    'petroRad_r',  # Petrosian radius in r-band (size)
                    'fracDeV_r',  # Fraction of light fit by de Vaucouleurs profile (bulge/disk indicator)
                    'expRad_r',  # Exponential fit radius (r-band)
                    'deVRad_r',  # de Vaucouleurs fit radius (r-band)
                    'expAB_r',  # Exponential fit axis ratio (r-band)
                    'deVAB_r',  # de Vaucouleurs fit axis ratio (r-band)
                    'expPhi_r',  # Exponential fit position angle (r-band)
                    'deVPhi_r',  # de Vaucouleurs fit position angle (r-band)
                    'psfMag_u', 'psfMag_g', 'psfMag_r', 'psfMag_i', 'psfMag_z',  # PSF magnitudes
                    'fiberMag_u', 'fiberMag_g', 'fiberMag_r', 'fiberMag_i', 'fiberMag_z',  # Fiber magnitudes
                    'cModelMag_u', 'cModelMag_g', 'cModelMag_r', 'cModelMag_i', 'cModelMag_z',  # Composite model magnitudes
                    'type',  # Object type (3 = galaxy)
                    'class',  # Spectral class (if available)
                    'clean',  # Clean photometry flag
                    'flags',  # Bitmask of photometric processing flags
                    'flags2',  # Additional photometric flags
                    'objc_type',  # Object classification (sometimes more detailed than type)
                    'rowc', 'colc',  # Row/column in SDSS image (for image cutouts)
                    'obj'  # Object number in field
                ]
            )
            
            if results is None or len(results) == 0:
                logger.warning("No galaxies found in the region")
                return []
            
            logger.info(f"Total objects found before filtering: {len(results)}")
            
            # Apply quality filters to ensure good data
            # Each filter is documented with its purpose
            mask = (
                (results['type'] == 3) &  # Galaxies only (type 3)
                (results['petroRad_r'] > 1) &  # Minimum size requirement
                (results['modelMagErr_g'] < 0.5) &  # Good g-band measurement
                (results['modelMagErr_r'] < 0.5)  # Good r-band measurement
            )
            
            # Log filtering statistics for debugging
            logger.info(f"Objects with type=3: {np.sum(results['type'] == 3)}")
            logger.info(f"Objects with clean=1: {np.sum(results['clean'] == 1)}")
            logger.info(f"Objects with petroRad_r>2: {np.sum(results['petroRad_r'] > 2)}")
            logger.info(f"Objects with modelMagErr_g<0.3: {np.sum(results['modelMagErr_g'] < 0.3)}")
            logger.info(f"Objects with modelMagErr_r<0.3: {np.sum(results['modelMagErr_r'] < 0.3)}")
            
            results = results[mask]
            logger.info(f"Total objects after filtering: {len(results)}")
            
            # Limit the number of results if specified
            if len(results) > limit:
                results = results[:limit]
            
            # Convert results to list of dictionaries with proper type handling
            galaxies = []
            for row in results:
                try:
                    # Calculate color indices for galaxy classification
                    g_r = (float(row['modelMag_g']) - float(row['modelMag_r'])) if row['modelMag_g'] is not None and row['modelMag_r'] is not None else None
                    u_r = (float(row['modelMag_u']) - float(row['modelMag_r'])) if row['modelMag_u'] is not None and row['modelMag_r'] is not None else None
                    
                    # Create standardized galaxy dictionary with all available data
                    galaxy = {
                        'objID': int(row['objID']),
                        'ra': float(row['ra']),
                        'dec': float(row['dec']),
                        'run': int(row['run']),
                        'rerun': int(row['rerun']),
                        'camcol': int(row['camcol']),
                        'field': int(row['field']),
                        'type': int(row['type']),
                        'class': str(row['class']) if row['class'] is not None else '',
                        'petroRad_r': float(row['petroRad_r']),
                        # Morphological parameters with null handling
                        'fracDeV_r': float(row['fracDeV_r']) if 'fracDeV_r' in row.colnames else None,
                        'expRad_r': float(row['expRad_r']) if 'expRad_r' in row.colnames else None,
                        'deVRad_r': float(row['deVRad_r']) if 'deVRad_r' in row.colnames else None,
                        'expAB_r': float(row['expAB_r']) if 'expAB_r' in row.colnames else None,
                        'deVAB_r': float(row['deVAB_r']) if 'deVAB_r' in row.colnames else None,
                        'expPhi_r': float(row['expPhi_r']) if 'expPhi_r' in row.colnames else None,
                        'deVPhi_r': float(row['deVPhi_r']) if 'deVPhi_r' in row.colnames else None,
                        # Model magnitudes with null handling
                        'modelMag': {
                            'u': float(row['modelMag_u']) if row['modelMag_u'] is not None else None,
                            'g': float(row['modelMag_g']) if row['modelMag_g'] is not None else None,
                            'r': float(row['modelMag_r']) if row['modelMag_r'] is not None else None,
                            'i': float(row['modelMag_i']) if row['modelMag_i'] is not None else None,
                            'z': float(row['modelMag_z']) if row['modelMag_z'] is not None else None
                        },
                        'modelMagErr': {
                            'u': float(row['modelMagErr_u']) if row['modelMagErr_u'] is not None else None,
                            'g': float(row['modelMagErr_g']) if row['modelMagErr_g'] is not None else None,
                            'r': float(row['modelMagErr_r']) if row['modelMagErr_r'] is not None else None,
                            'i': float(row['modelMagErr_i']) if row['modelMagErr_i'] is not None else None,
                            'z': float(row['modelMagErr_z']) if row['modelMagErr_z'] is not None else None
                        },
                        # PSF magnitudes with null handling
                        'psfMag': {
                            'u': float(row['psfMag_u']) if 'psfMag_u' in row.colnames and row['psfMag_u'] is not None else None,
                            'g': float(row['psfMag_g']) if 'psfMag_g' in row.colnames and row['psfMag_g'] is not None else None,
                            'r': float(row['psfMag_r']) if 'psfMag_r' in row.colnames and row['psfMag_r'] is not None else None,
                            'i': float(row['psfMag_i']) if 'psfMag_i' in row.colnames and row['psfMag_i'] is not None else None,
                            'z': float(row['psfMag_z']) if 'psfMag_z' in row.colnames and row['psfMag_z'] is not None else None
                        },
                        # Fiber magnitudes with null handling
                        'fiberMag': {
                            'u': float(row['fiberMag_u']) if 'fiberMag_u' in row.colnames and row['fiberMag_u'] is not None else None,
                            'g': float(row['fiberMag_g']) if 'fiberMag_g' in row.colnames and row['fiberMag_g'] is not None else None,
                            'r': float(row['fiberMag_r']) if 'fiberMag_r' in row.colnames and row['fiberMag_r'] is not None else None,
                            'i': float(row['fiberMag_i']) if 'fiberMag_i' in row.colnames and row['fiberMag_i'] is not None else None,
                            'z': float(row['fiberMag_z']) if 'fiberMag_z' in row.colnames and row['fiberMag_z'] is not None else None
                        },
                        # Composite model magnitudes with null handling
                        'cModelMag': {
                            'u': float(row['cModelMag_u']) if 'cModelMag_u' in row.colnames and row['cModelMag_u'] is not None else None,
                            'g': float(row['cModelMag_g']) if 'cModelMag_g' in row.colnames and row['cModelMag_g'] is not None else None,
                            'r': float(row['cModelMag_r']) if 'cModelMag_r' in row.colnames and row['cModelMag_r'] is not None else None,
                            'i': float(row['cModelMag_i']) if 'cModelMag_i' in row.colnames and row['cModelMag_i'] is not None else None,
                            'z': float(row['cModelMag_z']) if 'cModelMag_z' in row.colnames and row['cModelMag_z'] is not None else None
                        },
                        # Color indices
                        'g_r': g_r,
                        'u_r': u_r,
                        # Quality flags with null handling
                        'clean': bool(row['clean']) if 'clean' in row.colnames else None,
                        'flags': int(row['flags']) if 'flags' in row.colnames else None,
                        'flags2': int(row['flags2']) if 'flags2' in row.colnames else None,
                        'objc_type': int(row['objc_type']) if 'objc_type' in row.colnames else None,
                        # Image position with null handling
                        'rowc': float(row['rowc']) if 'rowc' in row.colnames else None,
                        'colc': float(row['colc']) if 'colc' in row.colnames else None,
                        'obj': int(row['obj']) if 'obj' in row.colnames else None
                    }
                    galaxies.append(galaxy)
                except (ValueError, TypeError) as e:
                    logger.error(f"Error processing galaxy data: {str(e)}")
                    continue
            
            logger.info(f"Found {len(galaxies)} galaxies matching the criteria")
            return galaxies
            
        except Exception as e:
            logger.error(f"Error querying galaxies: {str(e)}")
            return []

    def get_galaxy_images(self, galaxies: List[Dict]) -> List[str]:
        """
        Download color images for a list of galaxies using SDSS.
        Uses all available parameters from CSV to ensure precise image matching.
        Removes FITS files after conversion to JPG.
        
        Args:
            galaxies (List[Dict]): List of galaxy dictionaries containing ra and dec
            
        Returns:
            List[str]: List of paths to downloaded images
        """
        image_paths = []
        
        for galaxy in galaxies:
            try:
                # Calculate field size based on galaxy size
                field_size = min(max(galaxy['petroRad_r'] * 2, 0.3), 1.0)  # arcmin
                
                # Create coordinates object for precise positioning
                sky_coords = coords.SkyCoord(
                    ra=galaxy['ra'] * u.degree,
                    dec=galaxy['dec'] * u.degree,
                    frame='icrs'
                )
                
                # Create a precise query using all available parameters
                query = {
                    'run': galaxy['run'],
                    'rerun': galaxy['rerun'],
                    'camcol': galaxy['camcol'],
                    'field': galaxy['field'],
                    'data_release': 16,
                    'coordinates': sky_coords,
                    'radius': field_size * u.arcmin  # Use calculated field size
                }
                
                # Download g, r, i band images with precise parameters
                g_band = SDSS.get_images(
                    **query,
                    band='g'
                )
                
                r_band = SDSS.get_images(
                    **query,
                    band='r'
                )
                
                i_band = SDSS.get_images(
                    **query,
                    band='i'
                )
                
                # Verify images were downloaded
                if not g_band or not r_band or not i_band:
                    raise Exception("One or more images failed to download")
                
                # Define FITS file paths
                g_fits = f"data/images/{galaxy['objID']}_g.fits"
                r_fits = f"data/images/{galaxy['objID']}_r.fits"
                i_fits = f"data/images/{galaxy['objID']}_i.fits"
                
                try:
                    # Save raw FITS files with overwrite=True
                    g_band[0].writeto(g_fits, overwrite=True)
                    r_band[0].writeto(r_fits, overwrite=True)
                    i_band[0].writeto(i_fits, overwrite=True)
                    
                    # Get the image data from each band
                    # HDUList[0] is the primary HDU, and we need to access its data
                    g_data = g_band[0][0].data
                    r_data = r_band[0][0].data
                    i_data = i_band[0][0].data
                    
                    # Apply logarithmic scaling to enhance faint features
                    g_data = np.log10(g_data + 1)
                    r_data = np.log10(r_data + 1)
                    i_data = np.log10(i_data + 1)
                    
                    # Normalize each band using percentiles
                    def normalize_band(data):
                        data_min = np.percentile(data, 1)
                        data_max = np.percentile(data, 99)
                        return ((data - data_min) * (255.0 / (data_max - data_min))).clip(0, 255).astype(np.uint8)
                    
                    g_norm = normalize_band(g_data)
                    r_norm = normalize_band(r_data)
                    i_norm = normalize_band(i_data)
                    
                    # Create RGB image (using g, r, i bands for better color representation)
                    rgb_image = np.dstack((i_norm, r_norm, g_norm))  # i->R, r->G, g->B for better color balance
                    
                    # Create filename with galaxy properties
                    filename = f"galaxy_{galaxy['objID']}.jpg"
                    image_path = os.path.join(self.output_dir, "images", filename)
                    
                    # Create PIL Image and save as JPG
                    img = Image.fromarray(rgb_image)
                    img.save(image_path, quality=95)
                    
                    image_paths.append(image_path)
                    logger.info(f"Downloaded color image for galaxy {galaxy['objID']}")
                    
                finally:
                    # Clean up FITS files after processing
                    for fits_file in [g_fits, r_fits, i_fits]:
                        try:
                            if os.path.exists(fits_file):
                                os.remove(fits_file)
                                logger.debug(f"Removed temporary FITS file: {fits_file}")
                        except Exception as e:
                            logger.warning(f"Failed to remove FITS file {fits_file}: {str(e)}")
            
            except Exception as e:
                logger.error(f"Error downloading images for galaxy {galaxy['objID']}: {str(e)}")
                continue
        
        return image_paths

    def save_galaxy_metadata(self, galaxies, output_file):
        """
        Save galaxy metadata to a JSON file.
        
        Args:
            galaxies (List[Dict]): List of galaxy dictionaries
            output_file (str): Path to the output JSON file
        """
        with open(output_file, 'w') as f:
            json.dump(galaxies, f, indent=4)
        logger.info(f"Saved galaxy metadata to {output_file}")

    def process_galaxy_coordinates(self, csv_file: str, limit: int = 10) -> List[Dict]:
        """
        Read galaxy coordinates from a CSV file and download their images.
        Maps all relevant fields into the metadata dictionary for each galaxy.
        Args:
            csv_file (str): Path to the CSV file containing galaxy coordinates
            limit (int): Maximum number of galaxies to process. If None, process all.
        Returns:
            List[Dict]: List of galaxies with their properties
        """
        try:
            # Read the CSV file
            df = pd.read_csv(csv_file)
            logger.info(f"Read {len(df)} galaxies from {csv_file}")
            # Convert DataFrame to list of dictionaries
            galaxies = []
            for _, row in df.iterrows():
                try:
                    galaxy = {
                        'objID': int(row['objID']),
                        'ra': float(row['ra']),
                        'dec': float(row['dec']),
                        'run': int(row['run']),
                        'rerun': int(row['rerun']),
                        'camcol': int(row['camcol']),
                        'field': int(row['field']),
                        'modelMag_u': float(row['modelMag_u']) if not pd.isna(row.get('modelMag_u', None)) else None,
                        'modelMag_g': float(row['modelMag_g']) if not pd.isna(row.get('modelMag_g', None)) else None,
                        'modelMag_r': float(row['modelMag_r']) if not pd.isna(row.get('modelMag_r', None)) else None,
                        'modelMag_i': float(row['modelMag_i']) if not pd.isna(row.get('modelMag_i', None)) else None,
                        'modelMag_z': float(row['modelMag_z']) if not pd.isna(row.get('modelMag_z', None)) else None,
                        'modelMagErr_u': float(row['modelMagErr_u']) if not pd.isna(row.get('modelMagErr_u', None)) else None,
                        'modelMagErr_g': float(row['modelMagErr_g']) if not pd.isna(row.get('modelMagErr_g', None)) else None,
                        'modelMagErr_r': float(row['modelMagErr_r']) if not pd.isna(row.get('modelMagErr_r', None)) else None,
                        'modelMagErr_i': float(row['modelMagErr_i']) if not pd.isna(row.get('modelMagErr_i', None)) else None,
                        'modelMagErr_z': float(row['modelMagErr_z']) if not pd.isna(row.get('modelMagErr_z', None)) else None,
                        'psfMag_u': float(row['psfMag_u']) if not pd.isna(row.get('psfMag_u', None)) else None,
                        'psfMag_g': float(row['psfMag_g']) if not pd.isna(row.get('psfMag_g', None)) else None,
                        'psfMag_r': float(row['psfMag_r']) if not pd.isna(row.get('psfMag_r', None)) else None,
                        'psfMag_i': float(row['psfMag_i']) if not pd.isna(row.get('psfMag_i', None)) else None,
                        'psfMag_z': float(row['psfMag_z']) if not pd.isna(row.get('psfMag_z', None)) else None,
                        'fiberMag_u': float(row['fiberMag_u']) if not pd.isna(row.get('fiberMag_u', None)) else None,
                        'fiberMag_g': float(row['fiberMag_g']) if not pd.isna(row.get('fiberMag_g', None)) else None,
                        'fiberMag_r': float(row['fiberMag_r']) if not pd.isna(row.get('fiberMag_r', None)) else None,
                        'fiberMag_i': float(row['fiberMag_i']) if not pd.isna(row.get('fiberMag_i', None)) else None,
                        'fiberMag_z': float(row['fiberMag_z']) if not pd.isna(row.get('fiberMag_z', None)) else None,
                        'cModelMag_u': float(row['cModelMag_u']) if not pd.isna(row.get('cModelMag_u', None)) else None,
                        'cModelMag_g': float(row['cModelMag_g']) if not pd.isna(row.get('cModelMag_g', None)) else None,
                        'cModelMag_r': float(row['cModelMag_r']) if not pd.isna(row.get('cModelMag_r', None)) else None,
                        'cModelMag_i': float(row['cModelMag_i']) if not pd.isna(row.get('cModelMag_i', None)) else None,
                        'cModelMag_z': float(row['cModelMag_z']) if not pd.isna(row.get('cModelMag_z', None)) else None,
                        'petroRad_r': float(row['petroRad_r']) if not pd.isna(row.get('petroRad_r', None)) else None,
                        'fracDeV_r': float(row['fracDeV_r']) if not pd.isna(row.get('fracDeV_r', None)) else None,
                        'expRad_r': float(row['expRad_r']) if not pd.isna(row.get('expRad_r', None)) else None,
                        'deVRad_r': float(row['deVRad_r']) if not pd.isna(row.get('deVRad_r', None)) else None,
                        'expAB_r': float(row['expAB_r']) if not pd.isna(row.get('expAB_r', None)) else None,
                        'deVAB_r': float(row['deVAB_r']) if not pd.isna(row.get('deVAB_r', None)) else None,
                        'expPhi_r': float(row['expPhi_r']) if not pd.isna(row.get('expPhi_r', None)) else None,
                        'deVPhi_r': float(row['deVPhi_r']) if not pd.isna(row.get('deVPhi_r', None)) else None,
                        'rowc': float(row['rowc']) if not pd.isna(row.get('rowc', None)) else None,
                        'colc': float(row['colc']) if not pd.isna(row.get('colc', None)) else None,
                        'obj': int(row['obj']) if not pd.isna(row.get('obj', None)) else None,
                        'clean': int(row['clean']) if not pd.isna(row.get('clean', None)) else None,
                        'flags': int(row['flags']) if not pd.isna(row.get('flags', None)) else None,
                    }
                    galaxies.append(galaxy)
                except (ValueError, TypeError) as e:
                    logger.error(f"Error processing galaxy data: {str(e)}")
                    continue
            # Limit the number of results
            if limit is not None and len(galaxies) > limit:
                galaxies = galaxies[:limit]
            logger.info(f"Processing {len(galaxies)} galaxies from CSV")
            return galaxies
        except Exception as e:
            logger.error(f"Error reading galaxy coordinates from CSV: {str(e)}")
            return [] 

    def get_sdss_cutout_images(self, galaxies: List[Dict], scale: float = 0.396, width: int = 256, height: int = 256) -> List[str]:
        """
        Download SDSS cutout (postage stamp) images for a list of galaxies using the SDSS cutout service.
        Args:
            galaxies (List[Dict]): List of galaxy dictionaries containing ra, dec, objID
            scale (float): arcsec/pixel (default 0.396)
            width (int): width in pixels (default 64)
            height (int): height in pixels (default 64)
        Returns:
            List[str]: List of paths to downloaded cutout images
        """
        cutout_dir = os.path.join(self.output_dir, "sdss_cutouts")
        os.makedirs(cutout_dir, exist_ok=True)
        cutout_paths = []
        for galaxy in galaxies:
            ra = galaxy.get('ra')
            dec = galaxy.get('dec')
            objid = galaxy.get('objID')
            url = (
                f"https://skyserver.sdss.org/dr16/SkyServerWS/ImgCutout/getjpeg"
                f"?ra={ra}&dec={dec}&scale={scale}&width={width}&height={height}"
            )
            try:
                response = requests.get(url, timeout=10)
                if response.status_code == 200:
                    filename = os.path.join(cutout_dir, f"{objid}_cutout.jpg")
                    with open(filename, "wb") as f:
                        f.write(response.content)
                    cutout_paths.append(filename)
                    logger.info(f"Downloaded SDSS cutout for galaxy {objid}")
                else:
                    logger.warning(f"Failed to download cutout for {objid}: HTTP {response.status_code}")
            except Exception as e:
                logger.error(f"Error downloading SDSS cutout for {objid}: {str(e)}")
        return cutout_paths 

    def get_sdss_cutout_images_augmented(self, galaxies: List[Dict], n_versions: int = 3, shift_arcsec: float = 1.0, scale: float = 0.396, width: int = 256, height: int = 256) -> List[Dict]:
        """
        Download multiple SDSS cutouts per galaxy with small random center shifts, and save metadata for each cutout.
        Args:
            galaxies (List[Dict]): List of galaxy dictionaries containing ra, dec, objID
            n_versions (int): Number of cutouts per galaxy
            shift_arcsec (float): Max random shift in arcseconds
            scale (float): arcsec/pixel
            width (int): width in pixels
            height (int): height in pixels
        Returns:
            List[Dict]: List of metadata dicts for all cutouts (original fields + 'filename')
        """
        import random
        cutout_dir = os.path.join(self.output_dir, "sdss_cutouts")
        os.makedirs(cutout_dir, exist_ok=True)
        all_cutout_metadata = []
        for galaxy in galaxies:
            ra = galaxy.get('ra')
            dec = galaxy.get('dec')
            objid = galaxy.get('objID')
            for i in range(n_versions):
                # Random shift in arcseconds
                delta_ra = (random.uniform(-shift_arcsec, shift_arcsec) / 3600.0) / np.cos(np.deg2rad(dec))
                delta_dec = random.uniform(-shift_arcsec, shift_arcsec) / 3600.0
                shifted_ra = ra + delta_ra
                shifted_dec = dec + delta_dec
                url = (
                    f"https://skyserver.sdss.org/dr16/SkyServerWS/ImgCutout/getjpeg"
                    f"?ra={shifted_ra}&dec={shifted_dec}&scale={scale}&width={width}&height={height}"
                )
                try:
                    response = requests.get(url, timeout=10)
                    if response.status_code == 200:
                        filename = os.path.join(cutout_dir, f"{objid}_v{i}_cutout.jpg")
                        with open(filename, "wb") as f:
                            f.write(response.content)
                        logger.info(f"Downloaded SDSS cutout for galaxy {objid} version {i}")
                        # Save metadata for this cutout: only original fields + filename
                        cutout_metadata = dict(galaxy)
                        cutout_metadata['filename'] = filename
                        all_cutout_metadata.append(cutout_metadata)
                    else:
                        logger.warning(f"Failed to download cutout for {objid} version {i}: HTTP {response.status_code}")
                except Exception as e:
                    logger.error(f"Error downloading SDSS cutout for {objid} version {i}: {str(e)}")
        return all_cutout_metadata 

        """
        Create augmented versions of galaxy images for machine learning.
        
        This method:
        1. Creates multiple versions of each image with random transformations
        2. Applies consistent transformations across all bands
        3. Saves the augmented versions in separate directories
        
        Args:
            images (Dict[str, np.ndarray]): Dictionary of processed images for each band
            output_dir (str): Directory to save augmented versions
            num_versions (int): Number of augmented versions to create
            
        Note:
            The augmentation process:
            - Applies random rotations (0-360 degrees)
            - Applies random flips (horizontal and vertical)
            - Applies small random shifts
            - Maintains consistency across all bands
        """
        try:
            for i in range(num_versions):
                # Create version directory
                version_dir = os.path.join(output_dir, f"version_{i+1}")
                os.makedirs(version_dir, exist_ok=True)
                
                # Generate random transformations
                # Use same seed for all bands to maintain consistency
                np.random.seed(i)
                angle = np.random.uniform(0, 360)
                flip_h = np.random.choice([True, False])
                flip_v = np.random.choice([True, False])
                shift_x = np.random.randint(-5, 6)
                shift_y = np.random.randint(-5, 6)
                
                # Apply transformations to each band
                for band, image in images.items():
                    # Create augmented image
                    aug_image = image.copy()
                    
                    # Apply rotation
                    aug_image = rotate(aug_image, angle, reshape=False)
                    
                    # Apply flips
                    if flip_h:
                        aug_image = np.fliplr(aug_image)
                    if flip_v:
                        aug_image = np.flipud(aug_image)
                    
                    # Apply shift
                    aug_image = shift(aug_image, (shift_y, shift_x))
                    
                    # Save augmented image
                    output_path = os.path.join(version_dir, f"{band}.npy")
                    np.save(output_path, aug_image)
                
                logger.info(f"Created augmented version {i+1}")
            
        except Exception as e:
            logger.error(f"Error creating augmented versions: {str(e)}")