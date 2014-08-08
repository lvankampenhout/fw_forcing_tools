tools
=====

These scripts are from March 2014 and were used for realistic freshwater forcing experiments in CESM. 

* **find_bounding_cells.f90**: construct a mask of ocean cells around a target land mass. The mask is typically 1 cell thick.
* **def_blanking_mask.f90**: construct a mask of thickness > 1 around one or more landmasses. This mask can be used to blank existing values in the freshwater flux field before inserting your own (external) values. 
* Makefile
