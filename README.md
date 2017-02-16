tools
=====

The scripts in this folder are from March 2014 and helped doing our freshwater forcing experiments in CESM, doi:10.1002/2015GL064738

* **find_bounding_cells.f90**: construct a mask of ocean cells around a target land mass. These are the cells where the runoff is inserted. The mask is currently 1 cell wide. Note that this may be insufficient for high resolution experiments because of continuity restrictions, I have not explored this yet.
* **def_blanking_mask.f90**: construct a mask of thickness > 1 around one or more landmasses. This mask can be used to blank existing values in the freshwater flux field before inserting your own (external) values. 
* **Makefile**
