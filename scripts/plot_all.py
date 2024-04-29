"""
 Executes the following Python scripts to recreate figures from Savran et al.
  - pyCSEP: A Python Toolkit for Earthquake Forecast developers.
  Figure 1 and Figure 8 are not generated programatically, and are not included in the output from this script.

 This script can be executed in a stand-alone mode, but will require that the data files
 (10.5281/zenodo.5717419) are downloaded and extracted into the top-level directory of the
 reproducibility package. This script should be executed automatically when the Docker container
 containing the pyCSEP environment is created.

 Scripts executed :

     1. plot_Figure1.py
         Inputs:
             - No Input
         Outputs:
             - '../output/figure1.png'

     2. plot_Figure2.py
         Inputs:
             - '../data/area.dat'
             - '../data/GEAR1.dat'
             - '../data/NZHM_5year_rates-fromXML.dat'
  
         Outputs:
             - '../output/figure2.png'

     3. plot_Figure4-6.py
         Inputs:
             - '../data/area.dat'
             - '../data/GEAR1.dat'
             - '../data/NZHM_5year_rates-fromXML.dat'
             - '../data/nz5yrppe_c.dat'
             - '../data/nz5yrsup_c.dat'
             - '../data/GeoNet_catalog2021.txt'
         Outputs:
             - '../figures/figure4.png'
             - '../figures/figure4a.png'
             - '../figures/figure5.png'
             - '../figures/figure6a.png'
             - '../figures/figure6b.png'
             - '../figures/figure6c.png'

"""

import os
import time

import plot_fig_1
import plot_fig_2
import plot_fig_4to6


def verify_file_manifest():
    """ Checks directories for data and forecasts to determine which version of the reproducibility package to run.

        Returns:
            out (str): 'full' or 'light'
    """

    # Expected files for the 'full' and 'light' versions, full is light + full
    file_manifest = {
        'full': [
            '../forecasts/config.json',
            '../forecasts/m71_event.json',
            '../forecasts/results_complete.bin.gz'
        ]
    }

    print('Locating necessary files to recreate figures.')
    full = True
    for fpath in file_manifest['full']:
        if not os.path.exists(fpath):
            print(fpath)
            full = False
    # determine which version to run
    if full:
        output = 'full'
        print('Found all files, running full version of reproducibility package')
    else:
        raise FileNotFoundError('Missing files unable to run reproducibility package. '
                                'Try re-downloading from Zenodo. Contact k.graham [at] gns.cri.nz for assistance.')
    return output


def main(version):

    print(f'\n\nRunning {version} version of the reproducibility package. See README.md for more information.')
    print('=========================================================================================')

    print('')
    print('Generating Fig. 1')
    print('=================')
    plot_figure2.main()

    if ver == 'full':
        print('')
        print('Generating Fig. 2')
        print('=================')
        plot_fig_1.main()
    else:
        print("Skipping Fig. 3. See README for more information.")

    print('')
    print('Generating Fig. 4 to 6')
    print('=================')
    plot_fig_2.main()


    print('')
    print('Generating Fig. 5')
    print('=================')
    plot_fig_4to6.main()
    


if __name__ == "__main__":
    ver = verify_file_manifest()
    t0 = time.time()
    main(ver)
    t1 = time.time()
    print(f'Computed results in {t1 - t0:.3f} seconds.')
