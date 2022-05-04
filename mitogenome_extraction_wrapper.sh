#!/usr/bin/env bash
#
# Run Patchwork on the provided set of directories.

set -o errtrace
set -o errexit
set -o nounset
set -o pipefail

# Show line numbers and function names when debugging
export PS4='+ ${LINENO}:${FUNCNAME[0]:-}() '

readonly VERSION="0.1.0"
basedir=''
reference=''
files_only=0
id_offset=3
directory_list=''
kmer_offset=1
script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
main_script="${script_dir}/main.nf"

usage() { printf "%s" "\
usage:
  mitogenome_extraction_wrapper.sh [--help] [--version] DIRECTORY

description:
  Run the mitogenome extraction pipeline
  (github.com/ThiloSchulze/mitogenome-extraction) on the provided set of files.

recommended usage:
  - Display input files using '--files-only' to verify paths and species IDs
  - Once your happy with the input, do a normal run by removing the
    '--files-only' flag.
  - _If_ you need to perform multiple runs, set '--skip-combine'
  - After performing _multiple runs_, set '--combine-only' to combine outputs
    separately

options:
  miscellaneous:
    -h, --help          display this help message and exit
    -v, --version       display the version number and exit

  input/output control:
    -r, --reference     path to a FASTA file containing protein sequences to
                        match with (REQUIRED)

  offsets:
    -k, --kmer-offset     get K-mer size from the nth parent directory
                          (default: "$kmer_offset")
    -i, --id-offset       get ID from the nth parent directory (default: "$id_offset")

  flow control:
    -f, --files-only    display Patchwork input files and exit
"
  exit 1
}

# Display the version number and exit.
version_info() {
  echo "$VERSION"
  exit 1
}

# Display the provided error message before exiting with the provided status
# (default: 1).
error() {
  local message="$1"
  local status="${2-1}" # default exit status: 1
  echo "mitogenome_extraction_wrapper: error: $message"
  exit "$status"
}
export -f error

# Takes the path to a directory and an integer as an input. Go to the provided
# directory's parent directory X times.
up() {
  local dir="$1"
  local levels="$2"
  [[ ! -d "$dir" && ! -f "$dir" ]] && error "provided directory not found: $dir"
  while (( levels > 0 ))
  do
    dir=$( dirname "$dir" )
    (( levels = levels - 1 ))
  done
  echo "$dir"
}
export -f up

# Takes the path to a FASTQ file with forward reads and the path to another
# FASTQ file with reverse reads as an input. Returns a pattern that captures
# both of these reads. Also escapes dot to adhear to the expected Nextflow
# input format.
#
# example usage:
#   $ get_pattern NG-26745_1.fastq.gz NG-26745_2.fastq.gz
#   NG-26745_{1,2}\.fastq\.gz
get_pattern() {
  local file_a="$1"
  local file_b="$2"
  [[ "${#file_a}" -eq "${#file_b}" ]] || error "unequal length of filenames in \
raw reads: $file_a $file_b"
  mismatch_found=false
  before=''
  after=''
  pattern=''
  for (( i=0; i<${#file_a}; i++ ))
  do
    if [[ "${file_a:$i:1}" == "${file_b:$i:1}" ]]
    then
      if $mismatch_found
      then
        after+="${file_a:$i:1}"
      else
        before+="${file_a:$i:1}"
      fi
    else
      if $mismatch_found
      then
        error "found more than one difference in FASTQ files at position: $i"
      fi
      pattern="{${file_a:$i:1},${file_b:$i:1}}"
      mismatch_found=true
    fi
  done
  echo "${before}${pattern}${after}" | sed 's/\./\\\./g'
}
export -f get_pattern

# This function takes an array of assembly paths and returns the subset of
# those array with the smallest contigs.
get_smallest_contigs() {
  local assemblies=("${@}")
  local smallest_contigs=()
  local last_kmer_size=""
  local last_id=""
  local best_choice=""

  for assembly in "${assemblies[@]}"
  do
    kmer_size=$( basename "$( up "$assembly" "$kmer_offset" )" | tr -d 'K' )
    id=$( basename "$( up "$assembly" "$id_offset" )" )

    if [[ "$last_id" && "$last_kmer_size" ]]
    then
      if (( "$kmer_size" < "$last_kmer_size" ))
      then
        best_choice="$assembly"
      fi
      if [[ "$best_choice" && "$last_id" != "$id" ]]
      then
        smallest_contigs=("${smallest_contigs[@]}" "${best_choice}")
        best_choice="$assembly"
      fi
    else
      best_choice="$assembly"
    fi

    last_kmer_size="$kmer_size"
    last_id="$id"
  done
  if [[ "$best_choice" ]]
  then
    smallest_contigs=("${smallest_contigs[@]}" "${best_choice}")
  fi

  echo "${smallest_contigs[@]}"
}
export -f get_smallest_contigs

get_id() {
  local directory="${1}"
  local id_offset="${2}"
  basename "$( up "${directory}" "${id_offset}" )"
}
export -f get_id

get_outdir() {
  local directory="${1}"
  local id_offset="${2}"
  id="$( up "${directory}" "${id_offset}" )"
  dirout="${id}/mitogenome_extraction_out"
  echo "$dirout"
}
export -f get_outdir

get_raw_reads() {
  local directory="${1}"
  local id_offset="${2}"
  id="$( up "${directory}" "${id_offset}" )"
  raw_reads_path="${id}/raw_reads"
  [[ ! -d "$raw_reads_path" ]] && error "raw reads not found: $raw_reads_path"
  echo "$raw_reads_path"
}
export -f get_raw_reads

get_sequence_ids() {
  local filename="$1"
  local seq_ids=()
  readarray -d '' seq_ids < <( grep '^>' "$filename" | tr -d '>' )
  echo "${seq_ids[@]}"
}
export -f get_sequence_ids

run_mitogenome_extraction() {
  local species_id="$1"
  local contigs="$2"
  local reads_pattern="$3"
  local outdir="$4"
  local reference="$5"
  local main_script="$6"
  echo "$main_script\
 --mitogenome $reference\
 --contigs $contigs\
 --reads $reads_pattern\
 --species_id $species_id\
 --output $outdir\
 --max-cpus 4"
}
export -f run_mitogenome_extraction

run_on_contigs() {
  local reference="$1"
  local id_offset="$2"
  local files_only="$3"
  local main_script="$4"
  local contigs="$5"
  id="$( get_id "${contigs}" "${id_offset}" )"
  raw_reads_dir="$( get_raw_reads "${contigs}" "${id_offset}" )"
  reads_pattern="$( get_pattern $( fastq_files ${raw_reads_dir} 1 ) )"
  dirout="$( get_outdir "${contigs}" "${id_offset}" )"
  command="$( run_mitogenome_extraction "$id" "$contigs" "$reads_pattern"\
    "$dirout" "$reference" "$main_script" )"
  echo "ID:        $id"
  echo "Outdir:    ${dirout}"
  echo "Contigs:   ${contigs}"
  echo "Reference: ${reference}"
  echo "Reads:     ${reads_pattern}"
  echo -e "Command:   nextflow run ${command}\n"
  if (( ! files_only ))
  then
    mkdir -p "$dirout"
    # run command, write stdout + stderr to logfile
    ( cd "$dirout"; rm -f mitogenome_extraction.log;\
      nextflow run $command >> mitogenome_extraction.log 2>&1 )
  fi
}
export -f run_on_contigs

# Returns all FASTQ files found in the provided `path`
fastq_files() {
  local path="$1"
  find "$path" -type f -regex '\(.*fastq.gz\|.*fq.gz\)' | sort
  # find "$path" -type f \( --name "\*.fastq.gz" -o -name "\*.fq.gz" \) |\
  #   sort
}
export -f fastq_files

[[ $# -lt 1 ]] && usage

# Parse user-provided arguments
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
   -h | --help)
      usage
      ;;
    -v | --version)
      version_info
      ;;
    -r | --reference)
      reference="$2"
      shift # past argument
      shift # past value
      ;;
    -k | --kmer-offset)
      kmer_offset="$2"
      shift # past argument
      shift # past value
      ;;
    -f | --files-only)
      files_only=1
      shift # past argument
      ;;
    -i | --id-offset)
      id_offset="$2"
      shift # past argument
      shift # past value
      ;;
    -l | --from-list)
      directory_list="$2"
      shift # past argument
      shift # past value
      ;;
    --) # end of all options
      break
      ;;
    -*) # unknown option
      error "unknown option: $key"
      ;;
    *) # end of options
      [[ -n "$basedir" ]] && error "unrecognized positional argument: $key"
      basedir="$key"
      shift # past argument
      ;;
  esac
done

main() {
  # Display an error message if no arguments were provided
  [[ -z "$reference" ]] && error "missing mandatory argument --reference"
  reference="$( realpath "$reference" )"
  [[ ! -f "$reference" ]] && error "provided file not found: $reference"

  if [[ -n "$directory_list" ]]
  then # get assemblies from directory list
    smallest_contigs=()
    while read -r dir
    do
      readarray -d '' assemblies <\
        <(find "$dir" -type f -name "final_contigs.fasta" -not -wholename "*/work/*" -print0 )
      readarray smallest_contig<\
        <( get_smallest_contigs "${assemblies[@]}" )
      [[ -n "${smallest_contig[0]}" ]] && smallest_contigs+=("$smallest_contig")
    done < "$directory_list"
  else # find all assemblies in directory
    [[ -z "$basedir" ]] && error "missing mandatory argument DIRECTORY"
    [[ ! -d "$basedir" ]] && error "provided directory not found: $basedir"
    readarray -d '' assemblies <\
      <(find "$basedir" -type f -name "final_contigs.fasta" -not -wholename "*/work/*" -print0 )
    readarray smallest_contigs <\
      <( get_smallest_contigs "${assemblies[@]}" )
  fi

  if (( files_only ))
  then
    echo -e "These files will be used when \`--files-only\` is turned off:\n"
  fi

  for contigs in ${smallest_contigs[@]}
  do
    echo "Currently running: $contigs"
    run_on_contigs "$reference" "$id_offset" "$files_only"\
      "$main_script" "$contigs" || echo -e "Failed run: ${contigs}\n"; continue
  done

  echo 'Finished succesfully!'
  # To debug, run on the first file
  # run_on_contigs "$reference" "$id_offset" "$files_only"\
  #   "$main_script" ${smallest_contigs}
}

main
