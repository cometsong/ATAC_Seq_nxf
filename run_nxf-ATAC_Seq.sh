#!/bin/bash
##############################
  workflow="ATAC_Seq"
##############################
#TODO: implement generic nextflow workflow wrapper 'run' script; pass in workflow as $1 / build into run-WorkFlow scriptname?
#e.g. "run_NXF wkflw ARGS" => workflow="${1?-'Which Nextflow workflow would you like to run?'}"

##### CLI Args ###############
params=${1:-"$PWD/params.config"} && shift
profiles=${1:-"sumner"} && shift
args=${*} # e.g. '-resume'

##### Command Checks #########
has-cmd() { command -v "${1}" &>/dev/null; }
( has-cmd nextflow ) || module load nextflow \
  || { echo "Can't find command 'nextflow'" && exit 1; }
( has-cmd sbatch ) || module load slurm \
  || { echo "Can't find commands for 'slurm'" && exit 1; }

##### Get Datestamps #########
get_datestamp() { echo $(date +"${*:-"%Y-%m-%d %H:%M:%S"}"); } # from mbiome_utils.sh
today() { echo $(get_datestamp "%Y%m%d"); }

##### ENV Setup ##############
[ -n "${NXF_PIPES}" ] \
  || NXF_PIPES="/projects/${USER}/pipelines"

#TODO: ¿workflow=${this_script//*-} from name after hyphen?
#this_script=$(basename ${BASH_SOURCE[0]})

##### Workflow Path ##########
[ -d "${NXF_PIPES}" ] \
  && workflow_path="${NXF_PIPES}/${workflow}" \
  || workflow_path="${workflow}"

##### Run Name ###############
#NOTE: 'date +%6N' is nanoseconds =~ random-ish
run_name="${workflow}-$(date +%6N)"

##### Libdir #################
#NOTE: use below if libdir is centralized:
[ -d "${NXF_PIPES}/lib" ] \
  && libdir="${NXF_PIPES}/lib" \
  || libdir="${workflow_path}/lib"

##### tempdir Path ###########
# opt: tempdir="$PWD/${workflow,,}" # lowercase workflow name
# opt: tempdir=$(mktemp -d -p $TMPDIR -q nxf.XXXXX 2>/dev/null)||tempdir=$PWD
tempdir=$(mktemp -d -p /fastscratch -q nxf.XXXXX 2>/dev/null)||tempdir=$PWD

##### work-dir Path ##########
workdir="-work-dir ${tempdir}/NXF_work-${workflow}"
### Check args for work-dir ###
if [[ ${args[*]} =~ '-work-dir' ]] \
|| [[ ${args[*]} =~ '-w' ]] ; then
  workdir=''
fi

LOG="$PWD/nxf-$(today)-${run_name}.log"
#########################
##   The Main Action!  ##
nextflow                \
  -config ${params}     \
  run ${workflow_path}  \
  -profile ${profiles}  \
  -lib ${libdir}        \
  -name "${run_name}"   \
  ${workdir}            \
  ${args[*]}            \
2>&1 | tee -a ${LOG}
#########################
