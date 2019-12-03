# Converter from Gematronic vol data to wrfda input

Usage
=====

- Put all of your wrf_input in the WRF directory.
- Put all of your vol data in the VOL directory.
- In `vol2wrfdainput_gema.py`, change `os.environ['WRADLIB_DATA']` to your work directory.
- Change wrffile variable to the location of your wrf_input file (yes, step 1 is not mandatory, but let's be more tidy ;) )
- Run `vol2wrfdainput_gema.py` with command `python vol2wrfdainput_gema.py`.
- Your output will be in `out` directory with name `gemaob.radar`.
- Don't forget rename your `gemaob.radar` if you plan to re-run `vol2wrfdainput_gema.py` to process another domain.
