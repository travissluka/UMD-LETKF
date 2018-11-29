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

#determine the command used to run jobs (mpirun, aprun, neither...)
run_cmd=""
command -v mpirun && run_cmd="mpirun -n 2"
command -v aprun && run_cmd="aprun -n 2" 

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
$run_cmd $letkf_exe $json_file

# compare the results with the reference solution files
echo ""
echo "---------------------------------------------------"
echo " Checking output files against reference solutions"
echo "---------------------------------------------------"
echo ""

# check each file listed against the desired norm
for l in $(cat ref_solution/norms); do
    f=${l%%;*}
    n=${l##*;}
    e=${f##*.}

    echo "Checking difference in $f"
    
    if [[ ${e::4} == 'grib' ]]; then
	$check_grib_exe ref_solution/$f $f $n
    else
	$check_nc_exe ref_solution/$f $f $n
    fi
done
