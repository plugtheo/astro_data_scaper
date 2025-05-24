import json
import logging
import os

from galaxy_data_collector import GalaxyDataCollector

logger = logging.getLogger(__name__)

def main():
    # Initialize the collector
    collector = GalaxyDataCollector()
    
    # Process galaxies from CSV file
    # Get the root directory (one level up from src)
    root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    csv_file = os.path.join(root_dir, "data", "YOUR_CSV_FILE.csv")
    print(f"Looking for CSV file at: {csv_file}")  # Debug print
    galaxies = collector.process_galaxy_coordinates(csv_file, limit=None)
    
    if galaxies:
        # Download multiple SDSS cutout images per galaxy with augmentation
        cutout_metadata = collector.get_sdss_cutout_images_augmented(
            galaxies, n_versions=3, shift_arcsec=1.0, width=256, height=256
        )
        print(f"Downloaded {len(cutout_metadata)} SDSS cutout images (augmented)")
        
        # Save metadata for all cutouts
        cutout_metadata_path = os.path.join(collector.output_dir, "galaxy_cutout_metadata.json")
        collector.save_galaxy_metadata(cutout_metadata, cutout_metadata_path)
        print(f"Saved cutout metadata to {cutout_metadata_path}")
        
        # Save original metadata as well (optional, for reference)
        metadata_path = os.path.join(collector.output_dir, "galaxy_metadata.json")
        collector.save_galaxy_metadata(galaxies, metadata_path)
        print(f"Saved original metadata to {metadata_path}")
    else:
        print("No galaxies found")

if __name__ == "__main__":
    main() 