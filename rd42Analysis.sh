#!/bin/bash
# usage-message.sh

: ${1?"Usage: $0 source_files_dir settings_dir run_list run_number suffix"} ${2?"Usage: $0 source_files_dir settings_dir run_list run_number suffix"} ${3?"Usage: $0 source_files_dir settings_dir run_list run_number suffix"} ${4?"Usage: $0 source_files_dir settings_dir run_list run_number suffix"} ${5?"Usage: $0 source_files_dir settings_dir run_list run_number suffix"}

source_dir="$1"
settings_dir="$2"
out_dir="${source_dir}/output"
runlist="${3}"
run_no="$4"
suffix="$5"
scratch_out="/scratch/strip_telescope_tests/runDiego/output"

if [ ! -d "${source_dir}" ]; then
    echo "${source_dir} does not exist. Please give a valid source directory where the folders for the runs are"
    exit 1
fi

if [ ! -d "${source_dir}/${run_no}" ]; then
    echo "${source_dir} does not have run ${run_no}. Please give a valid run"
    exit 1
fi

if [ ! -d "${runlist}" ]; then
    echo "${runlist} does not exist. Please give the path to a valid run list directory"
    exit 1
fi

if [ ! -f "${runlist}/RunList_${run_no}.ini" ]; then
    echo "${runlist} does not have ${runlist}/RunList_${run_no}.ini"
    exit 1
fi

if [ ! -d "${settings_dir}" ]; then
    echo "${settings_dir} does not exist. Please give the path to a valid settings directory"
    exit 1
fi

if [ ! -f "${settings_dir}/settings.${run_no}.ini" ]; then
    echo "${settings_dir} does not have the settings file settings.${run_no}.ini for run ${run_no}"
    exit 1
fi

if [ ! -d "${out_dir}" ]; then
    echo "creating output directory ${out_dir}"
    mkdir "${out_dir}"
fi

if [ ! -d "${out_dir}/${suffix}" ]; then
    echo "creating output directory ${out_dir}/${suffix}"
    mkdir "${out_dir}/${suffix}"
fi

if [ ! -d "${out_dir}/${suffix}/${run_no}" ]; then
    mkdir "${out_dir}/${suffix}/${run_no}"
fi

if [ -L "${scratch_out}/${run_no}_${suffix}" ]; then
    rm "${scratch_out}/${run_no}_${suffix}"
fi
ln -s "${out_dir}/${suffix}/${run_no}" "${scratch_out}/${run_no}_${suffix}"


cp "${scratch_out}/index.php" "${out_dir}/${suffix}/${run_no}/"

diamondAnalysis -r "${runlist}/RunList_${run_no}.ini" -s "${settings_dir}" -o "${out_dir}/${suffix}" -i "${source_dir}/${run_no}"

exit 0

# ./rd42Analysis.sh /data/diamondSetup/diamondTestbeam/RD42-TestBeams/2016/cern_RD42_08_2016 ~/settings/full ~/RunLists 22022 full