#!/bin/bash
# PARALUS dependency installer
#
# Checks which PARALUS dependencies are available on this system and, after
# asking, installs all of them into a conda environment called "paralus".
# Supports Linux (x86_64) and Apple Silicon Macs.
#
# Usage: bash install_dependencies.sh [options]
#   --hisat2        install HISAT2 instead of STAR (low-RAM option)
#   --prefix PATH   create the environment at PATH instead of a named
#                   "paralus" environment (e.g. on HPC project storage)
#   --dry-run       show what would be done without installing anything
#   -h, --help      show this help

set -u

STAR_VERSION=2.7.11b
ENV_NAME=paralus
MIN_STAR_RAM_GB=32

use_hisat=0
prefix=""
dry_run=0

while [ $# -gt 0 ]; do
    case "$1" in
        --hisat2) use_hisat=1 ;;
        --prefix) shift; prefix="${1:-}"
                  [ -z "${prefix}" ] && { echo "error: --prefix needs a path"; exit 1; } ;;
        --dry-run) dry_run=1 ;;
        -h|--help) sed -n '2,14p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "error: unknown option $1 (see --help)"; exit 1 ;;
    esac
    shift
done

# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
run() {
    # runs a command, or only prints it with --dry-run
    if [ "${dry_run}" -eq 1 ]; then
        echo "  [dry-run] $*"
    else
        "$@"
    fi
}

ask() {
    # ask "question" default(y|n) -> returns 0 for yes
    local q="$1" def="$2" ans hint
    [ "${def}" = "y" ] && hint="[Y/n]" || hint="[y/N]"
    read -r -p "${q} ${hint} " ans
    ans="${ans:-${def}}"
    case "${ans}" in [Yy]*) return 0 ;; *) return 1 ;; esac
}

to_gb() { awk -v b="$1" 'BEGIN{printf "%d", b/1024/1024/1024}'; }

line() { echo "--------------------------------------------------------------------------"; }

# --------------------------------------------------------------------------
# 1. system
# --------------------------------------------------------------------------
os=$(uname -s); arch=$(uname -m)
if [ "${os}" = "Darwin" ] && [ "${arch}" = "arm64" ]; then
    system=mac
elif [ "${os}" = "Linux" ] && [ "${arch}" = "x86_64" ]; then
    system=linux
else
    echo "This installer supports Linux (x86_64) and Apple Silicon Macs."
    echo "Detected: ${os} ${arch}. Please follow INSTALLATION.md instead."
    exit 1
fi

is_hpc=0
if [ -n "${MODULESHOME:-}" ] || [ -n "${LMOD_CMD:-}" ] || command -v sbatch >/dev/null 2>&1; then
    is_hpc=1
fi

echo
echo "PARALUS dependency installer"
line
echo "System: $( [ "${system}" = mac ] && echo "Apple Silicon Mac" || echo "Linux x86_64" )"
[ "${dry_run}" -eq 1 ] && echo "Dry run: nothing will be installed."

if [ "${is_hpc}" -eq 1 ]; then
    echo
    echo "This looks like an HPC cluster. Clusters often provide the PARALUS"
    echo "dependencies as modules (e.g. module load star samtools bedtools2 cutadapt)."
    echo "If you can load all of them, you don't need this installer: just clone"
    echo "PARALUS (git clone https://github.com/KonstantinidouP/PARALUS.git) and"
    echo "load the modules before each run."
    echo
    ask "Continue and install the dependencies into a conda environment instead?" n || exit 0
fi

# --------------------------------------------------------------------------
# 2. what's already available
# --------------------------------------------------------------------------
echo
echo "Checking for PARALUS dependencies:"
missing=0
for tool in STAR hisat2 samtools bedtools cutadapt python3; do
    if p=$(command -v "${tool}" 2>/dev/null); then
        printf "  %-10s found   %s\n" "${tool}" "${p}"
    else
        printf "  %-10s missing\n" "${tool}"
        case "${tool}" in samtools|bedtools|cutadapt) missing=1 ;; esac
    fi
done
for mod in matplotlib numpy; do
    if python3 -c "import ${mod}" >/dev/null 2>&1; then
        printf "  %-10s found   (python module, optional - for PDF plots)\n" "${mod}"
    else
        printf "  %-10s missing (python module, optional - for PDF plots)\n" "${mod}"
    fi
done
if ! command -v STAR >/dev/null 2>&1 && ! command -v hisat2 >/dev/null 2>&1; then
    missing=1
fi

# --------------------------------------------------------------------------
# 3. memory and aligner choice
# --------------------------------------------------------------------------
if [ "${system}" = mac ]; then
    total_gb=$(to_gb "$(sysctl -n hw.memsize)")
    session_gb=${total_gb}
else
    total_gb=$(awk '/^MemTotal:/{printf "%d", $2/1024/1024}' /proc/meminfo)
    session_gb=${total_gb}
    # a job or session memory limit (SLURM / cgroups) can be lower than the machine's RAM
    for lim in /sys/fs/cgroup/memory.max /sys/fs/cgroup/memory/memory.limit_in_bytes; do
        if [ -r "${lim}" ]; then
            v=$(cat "${lim}")
            if [ "${v}" != "max" ]; then
                g=$(to_gb "${v}")
                [ "${g}" -lt "${session_gb}" ] && session_gb=${g}
            fi
        fi
    done
    if [ -n "${SLURM_MEM_PER_NODE:-}" ]; then
        g=$(( SLURM_MEM_PER_NODE / 1024 ))
        [ "${g}" -lt "${session_gb}" ] && session_gb=${g}
    fi
fi

echo
echo "Memory: ${total_gb} GB on this machine$( [ "${session_gb}" -lt "${total_gb}" ] && echo ", ${session_gb} GB available to this session" )"

if [ "${use_hisat}" -eq 0 ]; then
    if [ "${total_gb}" -lt "${MIN_STAR_RAM_GB}" ]; then
        echo
        echo "STAR needs at least ${MIN_STAR_RAM_GB} GB RAM for human or mouse genomes. On a computer"
        echo "with less RAM, PARALUS can use HISAT2 instead (run PARALUS with -H)."
        ask "Install HISAT2 instead of STAR (recommended for this computer)?" y && use_hisat=1
    elif [ "${session_gb}" -lt "${MIN_STAR_RAM_GB}" ]; then
        echo
        echo "Note: STAR needs at least ${MIN_STAR_RAM_GB} GB RAM for human or mouse genomes, and your"
        echo "current session has ${session_gb} GB. Request more memory for your session or job"
        echo "to run PARALUS with STAR. PARALUS can also work with HISAT2 (-H), but that's"
        echo "recommended only on computers with less than ${MIN_STAR_RAM_GB} GB RAM."
        echo "STAR will be installed. To install HISAT2 instead, rerun with --hisat2."
    fi
fi
aligner=$( [ "${use_hisat}" -eq 1 ] && echo hisat2 || echo STAR )
echo "Aligner to install: ${aligner}"

# --------------------------------------------------------------------------
# 4. what will be installed
# --------------------------------------------------------------------------
target=$( [ -n "${prefix}" ] && echo "${prefix}" || echo "conda environment '${ENV_NAME}'" )
echo
line
if [ "${missing}" -eq 0 ]; then
    echo "All required tools are already available, so PARALUS can run without a"
    echo "new environment. A separate environment is still useful to keep the"
    echo "PARALUS tools and versions together."
    ask "Create ${target} with all PARALUS dependencies anyway?" n || exit 0
else
    echo "The installer will create ${target} with all PARALUS"
    echo "dependencies: ${aligner}, samtools, bedtools, cutadapt, python, matplotlib, numpy."
    ask "Continue?" y || exit 0
fi

# Mac + STAR: reuse an existing compiled STAR, or compile it (needs Homebrew and gcc)
reuse_star=""
if [ "${system}" = mac ] && [ "${use_hisat}" -eq 0 ] && existing=$(command -v STAR 2>/dev/null); then
    ver=$(STAR --version 2>/dev/null)
    echo
    echo "Found STAR ${ver} at ${existing}"
    ask "Use this STAR (it will be copied into the new environment)? Answer n to compile a new one." y \
        && reuse_star="${existing}"
fi

if [ "${system}" = mac ] && [ "${use_hisat}" -eq 0 ] && [ -z "${reuse_star}" ]; then
    brew_bin=""
    for b in "$(command -v brew 2>/dev/null)" /opt/homebrew/bin/brew; do
        [ -n "${b}" ] && [ -x "${b}" ] && { brew_bin="${b}"; break; }
    done
    gxx=$(ls /opt/homebrew/bin/g++-[0-9]* 2>/dev/null | awk -F- '{print $NF" "$0}' | sort -n | tail -1 | cut -d" " -f2)
    echo
    echo "STAR has to be compiled on Apple Silicon Macs. This needs:"
    printf "  %-10s %s\n" "Homebrew" "$( [ -n "${brew_bin}" ] && echo "found" || echo "not found - will be installed (asks for your Mac password)" )"
    printf "  %-10s %s\n" "gcc" "$( [ -n "${gxx}" ] && echo "found (${gxx##*/})" || echo "not found - will be installed with Homebrew" )"
    echo "Compiling STAR takes a few minutes."
    if ! ask "Continue?" y; then
        if ask "Install HISAT2 instead of STAR?" y; then
            use_hisat=1; aligner=hisat2
        else
            echo "Stopped. See INSTALLATION.md to install STAR yourself."; exit 0
        fi
    fi
fi

# --------------------------------------------------------------------------
# 5. conda
# --------------------------------------------------------------------------
conda_bin=$(command -v conda 2>/dev/null)
if [ -z "${conda_bin}" ]; then
    for c in "$HOME/miniconda3/bin/conda" "$HOME/anaconda3/bin/conda"; do
        [ -x "${c}" ] && { conda_bin="${c}"; break; }
    done
fi
if [ -z "${conda_bin}" ]; then
    echo
    echo "conda was not found. It's needed to create the PARALUS environment."
    ask "Install Miniconda into $HOME/miniconda3 (no admin rights needed)?" y || { echo "Stopped."; exit 0; }
    if [ "${system}" = mac ]; then
        mc=Miniconda3-latest-MacOSX-arm64.sh
    else
        mc=Miniconda3-latest-Linux-x86_64.sh
    fi
    run curl -fL -o "/tmp/${mc}" "https://repo.anaconda.com/miniconda/${mc}" || { echo "error: download failed"; exit 1; }
    run bash "/tmp/${mc}" -b -p "$HOME/miniconda3" || { echo "error: Miniconda install failed"; exit 1; }
    run rm -f "/tmp/${mc}"
    run "$HOME/miniconda3/bin/conda" init "$(basename "${SHELL:-bash}")"
    conda_bin="$HOME/miniconda3/bin/conda"
fi

if [ -n "${prefix}" ]; then
    env_flag=(-p "${prefix}")
    env_path="${prefix}"
    activate="conda activate ${prefix}"
else
    env_flag=(-n "${ENV_NAME}")
    env_path="$("${conda_bin}" info --base 2>/dev/null)/envs/${ENV_NAME}"
    activate="conda activate ${ENV_NAME}"
fi

if [ -d "${env_path}" ]; then
    echo
    echo "An environment already exists at ${env_path}."
    ask "Install the PARALUS dependencies into it (existing packages are kept)?" n || { echo "Stopped."; exit 0; }
    conda_cmd=install
else
    conda_cmd=create
fi

pkgs="samtools bedtools cutadapt python matplotlib numpy"
[ "${use_hisat}" -eq 1 ] && pkgs="${pkgs} hisat2"

echo
echo "Installing with conda: ${pkgs} (this can take several minutes)"
# shellcheck disable=SC2086
run "${conda_bin}" "${conda_cmd}" -y "${env_flag[@]}" --override-channels \
    -c conda-forge -c bioconda ${pkgs} || { echo "error: conda ${conda_cmd} failed"; exit 1; }

# --------------------------------------------------------------------------
# 6. STAR
# --------------------------------------------------------------------------
if [ "${use_hisat}" -eq 0 ]; then
    build=$(mktemp -d)
    trap 'rm -rf "${build}"' EXIT
    if [ "${system}" = linux ]; then
        echo
        echo "Downloading the STAR ${STAR_VERSION} static binary"
        run curl -fL -o "${build}/star.zip" \
            "https://github.com/alexdobin/STAR/releases/download/${STAR_VERSION}/STAR_${STAR_VERSION}.zip" \
            || { echo "error: STAR download failed"; exit 1; }
        run "${env_path}/bin/python" -m zipfile -e "${build}/star.zip" "${build}"
        run cp "${build}/STAR_${STAR_VERSION}/Linux_x86_64_static/STAR" "${env_path}/bin/STAR"
        run chmod +x "${env_path}/bin/STAR"
    elif [ -n "${reuse_star}" ]; then
        echo
        echo "Copying your existing STAR into the environment"
        run cp "${reuse_star}" "${env_path}/bin/STAR"
    else
        if [ -z "${brew_bin}" ]; then
            echo
            echo "Installing Homebrew (it will ask for your Mac password)"
            run curl -fsSL -o "${build}/brew_install.sh" https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh \
                || { echo "error: Homebrew download failed"; exit 1; }
            run /bin/bash "${build}/brew_install.sh" || { echo "error: Homebrew install failed"; exit 1; }
            brew_bin=/opt/homebrew/bin/brew
        fi
        if [ -z "${gxx}" ]; then
            echo
            echo "Installing gcc with Homebrew"
            run "${brew_bin}" install gcc || { echo "error: gcc install failed"; exit 1; }
            gxx=$(ls /opt/homebrew/bin/g++-[0-9]* 2>/dev/null | awk -F- '{print $NF" "$0}' | sort -n | tail -1 | cut -d" " -f2)
            [ "${dry_run}" -eq 1 ] && [ -z "${gxx}" ] && gxx=/opt/homebrew/bin/g++-NN
        fi
        echo
        echo "Compiling STAR ${STAR_VERSION} (a few minutes; warnings are harmless)"
        run git clone --quiet https://github.com/alexdobin/STAR.git "${build}/STAR" \
            || { echo "error: could not download the STAR source"; exit 1; }
        run git -C "${build}/STAR" checkout --quiet "${STAR_VERSION}"
        run make -C "${build}/STAR/source" STARforMacStatic CXX="${gxx}" CXXFLAGS_SIMD= \
            >"${build}/make.log" 2>&1 \
            || { echo "error: STAR compilation failed. Last lines of the build log:"; tail -20 "${build}/make.log"; exit 1; }
        run cp "${build}/STAR/source/STAR" "${env_path}/bin/STAR"
    fi
fi

# --------------------------------------------------------------------------
# 7. check and finish
# --------------------------------------------------------------------------
echo
line
if [ "${dry_run}" -eq 1 ]; then
    echo "Dry run finished. Nothing was installed."
    exit 0
fi

echo "Checking the installed tools:"
ok=1
for tool in "${aligner}" samtools bedtools cutadapt; do
    if [ -x "${env_path}/bin/${tool}" ]; then
        printf "  %-10s ok\n" "${tool}"
    else
        printf "  %-10s NOT FOUND\n" "${tool}"; ok=0
    fi
done
if "${env_path}/bin/python" -c "import matplotlib, numpy" >/dev/null 2>&1; then
    printf "  %-10s ok\n" "python"
else
    printf "  %-10s plotting modules missing (PDF plots will be skipped)\n" "python"
fi

echo
if [ "${ok}" -eq 1 ]; then
    echo "PARALUS dependencies installed in ${env_path}"
    echo
    echo "Before each PARALUS run, activate the environment:"
    echo "  ${activate}"
    [ "${use_hisat}" -eq 1 ] && echo "and run PARALUS with the -H flag (HISAT2)."
else
    echo "Some tools are missing - see the messages above, or INSTALLATION.md."
    exit 1
fi
