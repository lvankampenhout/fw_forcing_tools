
FC=ifort
FFLAGS=-lnetcdf -lnetcdff
BIN=find_bounding_cells
BIN2=def_blanking_mask

all: find blank

find: 
	${FC} ${FFLAGS} find_bounding_cells.f90 -o ${BIN} -fpp

blank:
	${FC} ${FFLAGS} def_blanking_mask.f90 -o ${BIN2} -fpp


