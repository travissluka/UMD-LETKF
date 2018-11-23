#!/bin/bash
# This script launches the LETKF driver with a given json configuration file
# and then compares select output files with reference solutions to make sure
# they are within rounding errors of being identical

set -e

check_nc_exe=$1
check_grib_exe=$2
letkf_exe=$3
name=$4
config=$5
bin_dir=$6
src_dir=$7

# determine file locations
bkg_data=$bin_dir/$name.testdata
ref_data=$bin_dir/$name.ref_solutions/$config
wrk_dir=$bin_dir/$config
json_file=$config.json

# setup the working directory
if [ -d $wrk_dir ]; then
    rm $wrk_dir -r
fi
mkdir -p $wrk_dir
cd $wrk_dir
ln -sf $src_dir/* .
ln -sf $bkg_data/* .
mkdir -p ref_solution
ln -sf $ref_data/* ./ref_solution/

# run the LETKF
mpirun $letkf_exe $json_file

# compare the results with the reference solution files
echo ""
echo "---------------------------------------------------"
echo " Checking output files against reference solutions"
echo "---------------------------------------------------"
echo ""

shopt -s nullglob

# check the grib files
for f1 in ref_solution/*.{grib,grib2}; do
    f2=${f1##*/}
    echo "Checking difference in $f2"
    $check_grib_exe $f1 $f2 1e-14
done

# check the NetCDF files
for f1 in ref_solution/*.nc; do
    f2=${f1##*/}
    echo "Checking difference in $f2"
    $check_nc_exe $f1 $f2 1e-14
done
