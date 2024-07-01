#!/usr/bin/bash


########################################################################################################################
## Version: 1.1.1
## Author: Weining Li liwn@cau.edu.cn
## Date: 2023-07-05
##
## Calculate the average genetic distance for the provided populations
## 
## Usage: distance_multi_pop.sh --help
## 
## Required software/environment:
## R
## plink (1.90)
## 
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################


###################  Parameter processing  #####################
####################################################
## NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
## Parameters
TEMP=$(getopt -o h --long bfile:,out:,help \
              -n 'javawrap' -- "$@")
if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi
eval set -- "$TEMP"
## Parse parameters
while true; do
  case "$1" in
    --bfile )   bfile="$2";  shift 2 ;;  ## PLINK file prefixes, e.g.,"/public/home/popA /public/home/popB"
    --nchr )    nchr="$2" ;  shift 2 ;;  ## Number of chromosomes [30]
    --out )     out="$2" ;   shift 2 ;;  ## Output filename
    -h | --help )     grep " shift " $0 && exit 1 ;;
    -- ) shift; break ;;
    * ) shift; break ;;
  esac
done

## Default parameters
out=${out:="distance"}
nchr=${nchr:="30"}

## Check the scripts directory
if [[ ${code} ]]; then
  [[ ! -d ${code} ]] && echo "${code} not exists! " && exit 5
else
  script_path=$(dirname "$(readlink -f "$0")")
  code=$(dirname "$script_path")
fi

## Scripts
dist=${code}/R/mean_distance.R
func=${code}/shell/function.sh

## Load custom functions
[[ ! -s ${func} ]] && echo "Error: ${func} not found! " && exit 5
source ${func}

## Check if PLINK files exist
check_plink "${bfile}" ${nchr}

## Working directory
dir=$(dirname ${bfile})
cd ${dir} || exit

## Sort individuals by FID
plink --bfile ${bfile} --chr-set ${nchr} --indiv-sort natural --make-bed --out distance_tmp

## Sort individuals by FID
plink --bfile distance_tmp --chr-set ${nchr} --distance square0 1-ibs --out distance_tmp

## Calculate the mean of blocks corresponding to breed pairs in the matrix
$dist --prefix distance_tmp --out ${out}

## Remove intermediate files
rm distance_tmp.*
