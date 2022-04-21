COLIBREchemistry
=========

A python package that makes stellar abundance plots to analyse the COLIBRE simulations.

Requirements
----------------

The colibre-chemistry package requires:

+ `python3.6` or above
+ see requirements.txt

Usage
---------------

To run the script in the _single-run_ mode use
```bash
 python3 colibre_chemistry.py -d run_directory \
                              -s snapshot_name \
                              -c catalogue_name \
                              -n name_of_the_run \
                              -o path_to_output_directory 
```

To run the script in the _comparison_ mode use
```bash
 python3 colibre_chemistry.py -d directory_of_run1 directory_of_run2 \
                              -s snapshot_name_run1 snapshot_name_run2 \
                              -c catalogue_name_run1 catalogue_name_run2 \
                              -n name_of_the_run1 name_of_the_run2 \
                              -o path_to_output_directory
```



