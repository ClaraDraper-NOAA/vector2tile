From Mike, (Dec 1/2021): 

For now, I put the vector2tile tool here:

/scratch1/NCEPDEV/stmp2/Michael.Barlage/models/vector/vector2tile

I also put some test files here:

/scratch1/NCEPDEV/stmp2/Michael.Barlage/models/vector/v2t_data

Quick start from the vector2tile directory:

./configure
choose hera
load modules indicated
make to compile
./vector2tile_converter.exe namelist.vector2tile

the namelist defines the conversion direction and the paths of the files

the vector2tile pathway assumes that the vector file exists in the vector_restart_path directory and overwrites/creates tile files in the output_path

the tile2vector pathway is a little tricky, it assumes the tile files exist in tile_restart_path and overwrites only the snow variables in the vector file in the output_path, if the vector file does not exist in output_path the process will fail.

the overall assumption here is that we will have a full model vector restart file, then convert the vector to tiles for only snow variables, then convert the updated snow variables back to the full vector restart file

