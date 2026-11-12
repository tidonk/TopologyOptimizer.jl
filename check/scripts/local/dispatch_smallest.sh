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

# Usage: ./check/scripts/local/run_smallest.sh [--n <N>] [--settings <dir>] [--experiment_name <name>]


# ===== Path setup =====
script_dir=$(dirname "$(readlink -f "$0")")
project_dir="$(realpath "${script_dir}/../../..")"

INSTANCES_DIR="${project_dir}/check/instances/main"
SOLVE_SCRIPT="${project_dir}/check/scripts/cluster/solve_instance.jl"


# ===== Argument parsing =====
n_instances=10
settings_dir="${project_dir}/settings/subproblem_solving"
experiment_name="subproblem_solving_small"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --n)
            n_instances="$2"; shift 2 ;;
        --settings)
            settings_dir="$(realpath "$2")"; shift 2 ;;
        --experiment_name)
            experiment_name="$2"; shift 2 ;;
        *)
            echo "Unknown argument: $1" >&2
            echo "Usage: $0 [--n <N>] [--settings <dir>] [--experiment_name <name>]" >&2
            exit 1 ;;
    esac
done

# n smallest instances
mapfile -t casefiles < <(find "${INSTANCES_DIR}" -maxdepth 1 -type f -name "*.m" | sort -V | head -n "${n_instances}")
mapfile -t settings_files < <(find "${settings_dir}" -maxdepth 1 -name "*.toml" | sort)

echo "Instances (${#casefiles[@]}):"
printf '  %s\n' "${casefiles[@]##*/}"
echo "Settings (${#settings_files[@]}):"
printf '  %s\n' "${settings_files[@]##*/}"
echo "Results: ${project_dir}/check/results/${experiment_name}"


# ===== Main loop =====
for casefile in "${casefiles[@]}"; do
    for settings_file in "${settings_files[@]}"; do
        casename=$(basename "${casefile}" .m)
        settingsname=$(basename "${settings_file}" .toml)
        result_csv="${project_dir}/check/results/${experiment_name}/${settingsname}/${casename}.${settingsname}.csv"

        if [[ -f "${result_csv}" ]]; then
            echo "[Done   ] ${casename} / ${settingsname}"
            continue
        fi

        echo "[Solving] ${casename} / ${settingsname}"
        julia --startup-file=no --project="${project_dir}" "${SOLVE_SCRIPT}" "${casefile}" "${settings_file}" "${experiment_name}" \
            || echo "[Failed ] ${casename} / ${settingsname}" >&2
    done
done
