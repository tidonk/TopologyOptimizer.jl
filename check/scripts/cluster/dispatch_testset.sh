#!/bin/bash
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                           #
#  Copyright 2026 Operations Research, RWTH Aachen University               #
#                                                                           #
#  Licensed under the Apache License, Version 2.0 (the "License");          #
#  you may not use this file except in compliance with the License.         #
#  You may obtain a copy of the License at                                  #
#                                                                           #
#      http://www.apache.org/licenses/LICENSE-2.0                           #
#                                                                           #
#  Unless required by applicable law or agreed to in writing, software      #
#  distributed under the License is distributed on an "AS IS" BASIS,        #
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. #
#  See the License for the specific language governing permissions and      #
#  limitations under the License.                                           #
#                                                                           #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# ===== Path setup =====
script_dir=$(dirname "$(readlink -f "$0")")
project_dir="$(realpath "${script_dir}/../../..")"

INSTANCES_DIR="${project_dir}/check/instances/main"
RESULTS_DIR="${project_dir}/check/results"
LOG_FOLDER="${RESULTS_DIR}/cluster_logs"


# ===== Argument parsing =====
experiment_name="experiment"
SETTINGS_PATH=""
DRY_RUN=false
N_MIN_BUSES=73
N_MAX_BUSES=5000

while [[ $# -gt 0 ]]; do
    case "$1" in
        --experiment_name)
            experiment_name="$2"; shift 2 ;;
        --settings)
            SETTINGS_PATH="$2"; shift 2 ;;
        --dry-run)
            DRY_RUN=true; shift ;;
        --min_buses)
            N_MIN_BUSES="$2"; shift 2 ;;
        --max_buses)
            N_MAX_BUSES="$2"; shift 2 ;;
        *)
            echo "Unknown argument: $1" >&2
            echo "Usage: $0 --settings <file_or_dir> [--experiment_name <name>] [--min_buses <n>] [--max_buses <n>] [--dry-run]" >&2
            exit 1 ;;
    esac
done

if [[ -z "$SETTINGS_PATH" ]]; then
    echo "Error: --settings is required" >&2
    echo "Usage: $0 --settings <file_or_dir> [--experiment_name <name>] [--min_buses <n>] [--max_buses <n>] [--dry-run]" >&2
    exit 1
fi

SETTINGS_PATH="$(realpath "$SETTINGS_PATH")"

# Collect settings files: single .toml file or recursive directory scan
if [[ -f "$SETTINGS_PATH" ]]; then
    settings_files=("$SETTINGS_PATH")
    settings_desc="${SETTINGS_PATH} (1 file)"
elif [[ -d "$SETTINGS_PATH" ]]; then
    mapfile -t settings_files < <(find "${SETTINGS_PATH}" -name "*.toml" | sort)
    settings_desc="${SETTINGS_PATH} (${#settings_files[@]} files, recursive)"
else
    echo "Error: --settings must be an existing .toml file or directory: ${SETTINGS_PATH}" >&2
    exit 1
fi

if [[ ${#settings_files[@]} -eq 0 ]]; then
    echo "Error: no .toml files found at ${SETTINGS_PATH}" >&2
    exit 1
fi

dry_run_label=""
[[ "$DRY_RUN" == "true" ]] && dry_run_label=" (DRY RUN)"
echo "Dispatching experiment: ${experiment_name}${dry_run_label}"
echo "  Instances : ${INSTANCES_DIR}"
echo "  Settings  : ${settings_desc}"
echo "  Min buses : ${N_MIN_BUSES}"
echo "  Max buses : ${N_MAX_BUSES}"
echo "  Results   : ${RESULTS_DIR}/${experiment_name}"
echo "  Logs      : ${LOG_FOLDER}"

mkdir -p "${LOG_FOLDER}"
mkdir -p "${RESULTS_DIR}/temp"


# ===== Job counters =====
finished=0
running=0
submitted=0
resubmitted=0
unknown=0


# ===== Helper functions =====

# Print a colour-coded status line
# Usage: print_status <color> <status_label> <job_name> <message>
print_status() {
    local color="$1" status="$2" job_name="$3" message="$4"
    case "$color" in
        green)  echo -e "\e[32m[$(printf "%-10s" "$status") | $(printf "%-70s" "$job_name")] $message\e[0m" ;;
        yellow) echo -e "\e[33m[$(printf "%-10s" "$status") | $(printf "%-70s" "$job_name")] $message\e[0m" ;;
        blue)   echo -e "\e[34m[$(printf "%-10s" "$status") | $(printf "%-70s" "$job_name")] $message\e[0m" ;;
        red)    echo -e "\e[31m[$(printf "%-10s" "$status") | $(printf "%-70s" "$job_name")] $message\e[0m" ;;
    esac
}

# Return 0 (true) if a job with the given name is currently in the condor queue.
is_job_running() {
    echo "$_QUEUE_ARGS" | grep -qF "' '${1}' '"
}

# Refresh the condor queue cache
refresh_queue_cache() {
    _QUEUE_ARGS=$(condor_q -af Arguments 2>/dev/null)
}

# Return 0 (true) if the result CSV for this job already exists on disk
is_job_done() {
    [[ -f "$1" ]]
}

# Return 0 (true) if the job's stdout log exists and contains a Gurobi "Set parameter" line,
# indicating the solver was reached but the job terminated with unknown status (e.g. error/time-out)
has_unknown_termination() {
    [[ -f "$1" ]] && grep -q 'Set parameter' "$1"
}

# Extract hostname
get_last_execute_machine() {
    local log_file="$1"
    grep '^001 ' "$log_file" 2>/dev/null \
        | tail -n 1 \
        | sed 's/.*alias=\([^.&>]*\).*/\1/'
}

# Get the bus count from a case identifier (folder or basename)
# Supports names like "pglib_opf_case118_ieee" and "pglib_scopf_case67".
get_case_buses() {
    local identifier="$1"
    if [[ "$identifier" =~ _case([0-9]+) ]]; then
        echo "${BASH_REMATCH[1]}"
    else
        echo ""
    fi
}

# Write a HTCondor .sub file and submit it
# Usage: submit_condor_job <job_name> <sub_file_path> <casefile> <settings_file> <exp_name> <settings_name>
submit_condor_job() {
    local job_name="$1"
    local job_sub_path="$2"
    local casefile="$3"
    local settings_file="$4"
    local exp_name="$5"
    local settings_name="$6"

    cat <<SUBEOF > "${job_sub_path}"
# HTCondor job submission file for TopologyOptimizer experiment
universe = vanilla
executable = /bin/bash
arguments = "-c 'echo Running ${job_name} && mkdir -p /local/tmp/${USER} && $HOME/.juliaup/bin/julia --startup-file=no --history-file=no --project=${project_dir} ${script_dir}/solve_instance.jl ${casefile} ${settings_file} ${exp_name}'"
environment = "HOME=$HOME JULIA_PROJECT=${project_dir} GRB_LICENSE_FILE=${GRB_LICENSE_FILE:?Set GRB_LICENSE_FILE before submitting jobs} JULIA_DEPOT_PATH=/local/tmp/${USER}/.julia"

# Log and output file paths
output = ${LOG_FOLDER}/${job_name}.out
error  = ${LOG_FOLDER}/${job_name}.err
log    = ${LOG_FOLDER}/${job_name}.log

# Job details
batch_name = ${exp_name}: ${settings_name}
name = ${job_name}

# Automatically remove jobs running longer than 4 hours (14400s)
periodic_remove = (JobStatus == 2 && (CurrentTime - JobStartDate) > 14400)

# Exclude the faulty machine
queue
SUBEOF

    [[ "$DRY_RUN" == "false" ]] && condor_submit "${job_sub_path}" > /dev/null 2>&1
}


# ===== Cache condor queue state & show pre-flight summary =====
refresh_queue_cache
queue_total=$(echo "$_QUEUE_ARGS" | grep -c "Running ${experiment_name}_" 2>/dev/null || true)
echo ""
echo "Condor queue: ${queue_total} job(s) for experiment '${experiment_name}' already queued"
echo ""

# ===== Main dispatch loop =====
echo "Submitting jobs..."

mapfile -t casefolder_paths < <(find "${INSTANCES_DIR}" -maxdepth 1 -mindepth 1 -type d | sort -V)
mapfile -t direct_casefiles < <(find "${INSTANCES_DIR}" -maxdepth 1 -mindepth 1 -type f -name "*.m" | sort -V)

process_casefile() {
    local casefile="$1"
    local casefolder="$2" # for batch_name; can be "local" for root-level files

    case_basename=$(basename "$casefile")
    casename="${case_basename%.m}"

    n_buses=$(get_case_buses "$casename")
    if [[ -n "$n_buses" && "$n_buses" -lt "$N_MIN_BUSES" ]]; then
        echo "  [Skipping ${casefile}: ${n_buses} buses < ${N_MIN_BUSES}]"
        return
    fi
    if [[ -n "$n_buses" && "$n_buses" -gt "$N_MAX_BUSES" ]]; then
        echo "  [Skipping ${casefile}: ${n_buses} buses > ${N_MAX_BUSES}]"
        return
    fi

    for settings_file in "${settings_files[@]}"; do
        settings_basename=$(basename "$settings_file")
        settingsname="${settings_basename%.toml}"

        JOB_NAME="${experiment_name}_${casename}.${settingsname}"
        JOB_SUB_PATH="${RESULTS_DIR}/temp/${JOB_NAME}.sub"
        RESULT_CSV="${RESULTS_DIR}/${experiment_name}/${settingsname}/${casename}.${settingsname}.csv"
        OUT_FILE="${LOG_FOLDER}/${JOB_NAME}.out"

        if is_job_running "$JOB_NAME"; then
            print_status "yellow" "Running" "$JOB_NAME" "Job is already running"
            echo "            | Log: ${LOG_FOLDER}/${JOB_NAME}.log"
            running=$((running + 1))
            continue
        fi

        if is_job_done "$RESULT_CSV"; then
            print_status "green" "Done" "$JOB_NAME" "Result CSV already exists"
            if is_job_running "$JOB_NAME"; then
                print_status "red" "Redundant" "$JOB_NAME" "Done but still running in queue"
                running=$((running + 1))
            fi
            finished=$((finished + 1))
            continue
        fi

        PRIOR_ERROR=false
        if has_unknown_termination "$OUT_FILE"; then
            PRIOR_ERROR=true
            unknown=$((unknown + 1))
        fi

        if [[ "$PRIOR_ERROR" == "true" ]]; then
            last_machine=$(get_last_execute_machine "${LOG_FOLDER}/${JOB_NAME}.log")
            if [[ "$DRY_RUN" == "true" ]]; then
                print_status "red" "Dry Resub" "$JOB_NAME" "Would resubmit (prior error on ${last_machine:-unknown})"
                echo "            | Sub: ${JOB_SUB_PATH}"
            else
                print_status "red" "Resubmit" "$JOB_NAME" "Resubmitting after prior error on ${last_machine:-unknown}"
                echo "            | Err: ${LOG_FOLDER}/${JOB_NAME}.err"
                echo "            | Out: ${OUT_FILE}"
            fi
            resubmitted=$((resubmitted + 1))
        else
            if [[ "$DRY_RUN" == "true" ]]; then
                print_status "yellow" "Dry Run" "$JOB_NAME" "Would submit — .sub file written for inspection"
                echo "            | Sub: ${JOB_SUB_PATH}"
            else
                print_status "blue" "Submitting" "$JOB_NAME" "Submitting condor job"
                echo "            | Log: ${LOG_FOLDER}/${JOB_NAME}.log"
            fi
            submitted=$((submitted + 1))
        fi

        submit_condor_job \
            "$JOB_NAME" "$JOB_SUB_PATH" \
            "$casefile" "$settings_file" \
            "$experiment_name" "$settingsname"
    done
}

for casefolder_path in "${casefolder_paths[@]}"; do
    casefolder=$(basename "$casefolder_path")

    # Extract bus count from folder name (e.g. pglib_opf_case118_ieee -> 118)
    n_buses=$(get_case_buses "$casefolder")
    if [[ -n "$n_buses" && "$n_buses" -lt "$N_MIN_BUSES" ]]; then
        echo "  [Skipping casefolder ${casefolder}: ${n_buses} buses < ${N_MIN_BUSES}]"
        continue
    fi
    if [[ -n "$n_buses" && "$n_buses" -gt "$N_MAX_BUSES" ]]; then
        echo "  [Skipping casefolder ${casefolder}: ${n_buses} buses > ${N_MAX_BUSES}]"
        continue
    fi

    mapfile -t casefiles < <(find "${casefolder_path}" -maxdepth 1 -name "*.m" | sort -V)

    for casefile in "${casefiles[@]}"; do
        process_casefile "$casefile" "$casefolder"
    done
done

for casefile in "${direct_casefiles[@]}"; do
    process_casefile "$casefile" "root"
done


# ===== Summary =====
printf "\n\e[4mSummary:\e[0m\n"
printf "\e[0m%-10s \e[32m%-10s \e[33m%-10s \e[34m%-10s \e[31m%-10s \e[31m%-10s\e[0m\n" \
    "Status" "Finished" "Running" "Submitted" "Resubmit" "Unknown"
printf "%-10s \e[32m%-10s \e[33m%-10s \e[34m%-10s \e[31m%-10s \e[31m%-10s\e[0m\n" \
    " " "$finished" "$running" "$submitted" "$resubmitted" "$unknown"
