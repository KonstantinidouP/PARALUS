#!/bin/bash
# PARALUS installer
#
# Checks which PARALUS dependencies are available on this system and, after
# asking, installs all of them into a conda environment called "paralus".
# It also sets up the "paralus" command inside that environment.
# Supports Linux (x86_64) and Apple Silicon Macs.
#
# Usage: bash install_paralus.sh [options]
#   --hisat2          install HISAT2 instead of STAR (low-RAM option)
#   --env NAME        use the conda environment NAME instead of "paralus",
#                     e.g. an environment you already have
#   --prefix PATH     create the environment at PATH instead of a named
#                     "paralus" environment (e.g. on HPC project storage)
#   --launcher-only   only set up the "paralus" command (adds the PARALUS
#                     folder to your PATH), e.g. if you load the
#                     dependencies as HPC modules
#   --dry-run         show what would be done without installing anything
#   -h, --help        show this help

set -u

STAR_VERSION=2.7.11b
ENV_NAME=paralus
MIN_STAR_RAM_GB=32

use_hisat=0
prefix=""
env_given=0
dry_run=0
launcher_only=0

while [ $# -gt 0 ]; do
    case "$1" in
        --hisat2) use_hisat=1 ;;
        --env) shift; ENV_NAME="${1:-}"; env_given=1
               [ -z "${ENV_NAME}" ] && { echo "error: --env needs an environment name"; exit 1; } ;;
        --prefix) shift; prefix="${1:-}"
                  [ -z "${prefix}" ] && { echo "error: --prefix needs a path"; exit 1; } ;;
        --launcher-only) launcher_only=1 ;;
        --dry-run) dry_run=1 ;;
        -h|--help) sed -n '2,19p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "error: unknown option $1 (see --help)"; exit 1 ;;
    esac
    shift
done
if [ "${env_given}" -eq 1 ] && [ -n "${prefix}" ]; then
    echo "error: use either --env or --prefix, not both"; exit 1
fi

# the PARALUS folder is the folder this installer is in
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ ! -f "${here}/paralus" ]; then
    echo "error: the paralus launcher was not found in ${here}"
    echo "Run this installer from inside the PARALUS folder, where it came from."
    exit 1
fi

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

next_steps_footer() {
    echo
    echo "If you move the PARALUS folder, run this installer again from its new location."
}

MARKER="# Added for PARALUS - delete this line and the next one to undo"

shell_rc() {
    # the startup file of the user's shell
    case "$(basename "${SHELL:-bash}")" in
        zsh) echo "$HOME/.zshrc" ;;
        *)   echo "$HOME/.bashrc" ;;
    esac
}

tilde() {
    # shows a path in the home folder as ~/... (works with the bash 3.2 of macOS)
    case "$1" in
        "$HOME"/*) echo "~${1#"$HOME"}" ;;
        *)         echo "$1" ;;
    esac
}

check_login_file() {
    # On ssh logins (e.g. HPC clusters) bash reads ~/.bash_profile (or
    # ~/.bash_login or ~/.profile) instead of ~/.bashrc. Warn if that file
    # does not load ~/.bashrc, since the PATH line would then be skipped.
    local f login=""
    for f in "$HOME/.bash_profile" "$HOME/.bash_login" "$HOME/.profile"; do
        [ -f "${f}" ] && { login="${f}"; break; }
    done
    if [ -n "${login}" ] && grep -q "bashrc" "${login}"; then
        return 0
    fi
    echo
    if [ -n "${login}" ]; then
        echo "Note: when you log in remotely (e.g. over ssh), bash reads $(tilde "${login}")"
        echo "instead of ~/.bashrc, and $(tilde "${login}") does not load ~/.bashrc. To have"
        echo "paralus available after every login, add this line to $(tilde "${login}"):"
    else
        echo "Note: when you log in remotely (e.g. over ssh), bash reads ~/.bash_profile"
        echo "instead of ~/.bashrc, and you don't have a ~/.bash_profile. To have"
        echo "paralus available after every login, create ~/.bash_profile with this line:"
    fi
    echo '  [ -f ~/.bashrc ] && . ~/.bashrc'
}

old_paralus_lines() {
    # prints the PARALUS marker line(s) in a startup file, each with the
    # line after it only if that line is a PATH export
    awk -v m="${MARKER}" '
        pending && /^export PATH=".*:\$PATH"$/ { print; pending=0; next }
        { pending=0 }
        $0==m { print; pending=1 }' "$1"
}

remove_paralus_lines() {
    # removes the PARALUS marker line(s) from a startup file. The line after a
    # marker is removed only if it is a PATH export, so nothing else is lost
    # if that line was edited or deleted. A backup is kept as <file>.paralus_backup
    local rc="$1"
    cp "${rc}" "${rc}.paralus_backup" || return 1
    awk -v m="${MARKER}" '
        pending && /^export PATH=".*:\$PATH"$/ { pending=0; next }
        { pending=0 }
        $0==m { pending=1; next }
        { print }' "${rc}.paralus_backup" > "${rc}.paralus_tmp" \
        && cat "${rc}.paralus_tmp" > "${rc}" && rm -f "${rc}.paralus_tmp"
}

setup_path() {
    # adds the PARALUS folder to PATH in the user's shell startup file
    # (pass "noask" when the user has already agreed)
    local confirm="${1:-ask}" rc
    rc=$(shell_rc)
    local newline="export PATH=\"${here}:\$PATH\""

    echo
    if [ -f "${rc}" ] && grep -qxF "${newline}" "${rc}"; then
        echo "PARALUS is already on your PATH in ${rc}"
    else
        echo "This adds the PARALUS folder to your PATH in ${rc}, so you can run"
        echo "paralus from any folder. Nothing else will be installed."
        if [ "${confirm}" = ask ]; then
            ask "Continue?" y || { echo "Stopped."; exit 0; }
        fi
        if [ "${dry_run}" -eq 1 ]; then
            echo "  [dry-run] add to ${rc}: ${newline}"
        else
            if [ -f "${rc}" ] && grep -qxF "${MARKER}" "${rc}"; then
                # an earlier PARALUS line (e.g. from a folder that was moved): replace it
                remove_paralus_lines "${rc}"
            fi
            printf '%s\n%s\n' "${MARKER}" "${newline}" >> "${rc}"
        fi
        echo
        echo "Added to ${rc}:"
        echo "  ${newline}"
    fi
    run chmod +x "${here}/paralus"
    [ "${rc}" = "$HOME/.bashrc" ] && check_login_file

    echo
    line
    echo "Next steps:"
    echo "  1. Reload your shell settings:  source $(tilde "${rc}")"
    echo "  2. Make sure the PARALUS dependencies are available"
    echo "     (e.g. load their modules on an HPC cluster)."
    echo "  3. Check that it works:        paralus --help"
    next_steps_footer
}

# --------------------------------------------------------------------------
# --launcher-only
# --------------------------------------------------------------------------
if [ "${launcher_only}" -eq 1 ]; then
    echo
    echo "PARALUS launcher setup"
    line
    [ "${dry_run}" -eq 1 ] && echo "Dry run: nothing will be changed."
    setup_path
    exit 0
fi

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
echo "PARALUS installer"
line
echo "System: $( [ "${system}" = mac ] && echo "Apple Silicon Mac" || echo "Linux x86_64" )"
[ "${dry_run}" -eq 1 ] && echo "Dry run: nothing will be installed."

# --------------------------------------------------------------------------
# conda and the target environment (conda itself is installed later if needed)
# --------------------------------------------------------------------------
conda_bin=$(command -v conda 2>/dev/null)
if [ -z "${conda_bin}" ]; then
    for c in "$HOME/miniconda3/bin/conda" "$HOME/anaconda3/bin/conda"; do
        [ -x "${c}" ] && { conda_bin="${c}"; break; }
    done
fi

set_env_path() {
    if [ -n "${prefix}" ]; then
        env_flag=(-p "${prefix}")
        env_path="${prefix}"
        env_desc="${prefix}"
        activate="conda activate ${prefix}"
    else
        env_flag=(-n "${ENV_NAME}")
        # where conda keeps this environment (base, ~/.conda/envs on clusters, ...);
        # a new environment goes to the usual envs folder
        env_path=$("${conda_bin}" env list 2>/dev/null \
            | awk -v n="${ENV_NAME}" '$1==n {print $NF; exit}')
        [ -n "${env_path}" ] || env_path="$("${conda_bin}" info --base 2>/dev/null)/envs/${ENV_NAME}"
        env_desc="the conda environment '${ENV_NAME}'"$'\n'"(${env_path})"
        activate="conda activate ${ENV_NAME}"
    fi
}
env_path=""
[ -n "${conda_bin}" ] || [ -n "${prefix}" ] && set_env_path

link_launcher() {
    echo
    echo "Setting up the paralus command"
    run chmod +x "${here}/paralus"
    run ln -sf "${here}/paralus" "${env_path}/bin/paralus"

    # PARALUS folder added to PATH by hand (older README) or with --launcher-only
    local rc found
    for rc in "$HOME/.zshrc" "$HOME/.bashrc"; do
        [ -f "${rc}" ] || continue
        found=$(old_paralus_lines "${rc}")
        [ -n "${found}" ] || continue
        echo
        echo "Found an earlier PARALUS setup in ${rc}:"
        echo "${found}" | sed 's/^/  /'
        echo "It may point to an older PARALUS version. The paralus command now comes"
        echo "with the conda environment, so these lines are no longer needed."
        if ask "Remove them?" y; then
            if [ "${dry_run}" -eq 1 ]; then
                echo "  [dry-run] remove these lines from ${rc}"
            elif remove_paralus_lines "${rc}"; then
                echo "Removed. The previous version of the file is saved as ${rc}.paralus_backup"
            else
                echo "Could not change ${rc}. You can delete these lines yourself."
            fi
        fi
    done
}

finish() {
    # final check and message for a conda environment install
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
    if [ -x "${env_path}/bin/paralus" ]; then
        printf "  %-10s ok\n" "paralus"
    else
        printf "  %-10s NOT FOUND\n" "paralus"; ok=0
    fi

    echo
    if [ "${ok}" -eq 1 ]; then
        echo "PARALUS dependencies installed in ${env_desc}"
        echo
        echo "Next steps:"
        echo "  1. Reload your shell settings:  source $(tilde "$(shell_rc)")"
        echo "     (or close this terminal window and open a new one)"
        printf "  2. Activate the environment:   %s\n" "${activate}"
        echo "  3. Check that it works:        paralus --help"
        echo
        echo "Activate the environment before each PARALUS run."
        [ "${aligner}" = hisat2 ] && echo "Run PARALUS with the -H flag, since HISAT2 was installed instead of STAR."
        if [ -n "${prefix}" ]; then
            rerun="bash install_paralus.sh --prefix ${prefix}"
        elif [ "${env_given}" -eq 1 ]; then
            rerun="bash install_paralus.sh --env ${ENV_NAME}"
        else
            rerun=""
        fi
        if [ -n "${rerun}" ]; then
            echo "If you move the PARALUS folder, run this installer again from its new location:"
            echo "  ${rerun}"
        else
            echo "If you move the PARALUS folder, run this installer again from its new location."
        fi
        exit 0
    else
        echo "Some tools are missing - see the messages above, or INSTALLATION.md."
        exit 1
    fi
}

# An environment that already has everything (e.g. rerun after moving the
# PARALUS folder): only update the paralus command.
shortcut_if_complete() {
    [ -n "${env_path}" ] && [ -d "${env_path}" ] || return 0
    local has=1 tool
    for tool in samtools bedtools cutadapt; do
        [ -x "${env_path}/bin/${tool}" ] || has=0
    done
    if [ -x "${env_path}/bin/STAR" ] && [ "${use_hisat}" -eq 0 ]; then
        aligner=STAR
    elif [ -x "${env_path}/bin/hisat2" ]; then
        aligner=hisat2
    else
        has=0
    fi
    if [ "${has}" -eq 1 ]; then
        echo
        echo "Found an existing environment with all PARALUS dependencies:"
        echo "  ${env_path}"
        echo "Only the paralus command will be updated."
        link_launcher
        finish
    fi
}
shortcut_if_complete

if [ "${is_hpc}" -eq 1 ]; then
    echo
    echo "This looks like an HPC cluster. Clusters often provide the PARALUS"
    echo "dependencies as modules (e.g. module load star samtools bedtools2 cutadapt)."
    echo "If you can load all of them, you don't need conda. Run this instead to set"
    echo "up the paralus command only:"
    echo "  bash install_paralus.sh --launcher-only"
    echo "and load the modules before each run."
    echo
    ask "Continue and install the dependencies into a conda environment instead?" n || exit 0
fi

# --------------------------------------------------------------------------
# 2. what's already available
# --------------------------------------------------------------------------
# An existing target environment is checked itself; otherwise the tools
# already available in this shell are checked.
env_exists=0
[ -n "${env_path}" ] && [ -d "${env_path}" ] && env_exists=1
find_tool() {
    if [ "${env_exists}" -eq 1 ]; then
        [ -x "${env_path}/bin/$1" ] && echo "${env_path}/bin/$1"
    else
        command -v "$1" 2>/dev/null
    fi
}
if [ "${env_exists}" -eq 1 ]; then
    py="${env_path}/bin/python3"
else
    py=python3
fi

echo
if [ "${env_exists}" -eq 1 ]; then
    echo "Checking for PARALUS dependencies in ${env_path}:"
else
    echo "Checking for PARALUS dependencies:"
fi
missing=0
for tool in STAR hisat2 samtools bedtools cutadapt python3; do
    if p=$(find_tool "${tool}"); then
        printf "  %-10s found   %s\n" "${tool}" "${p}"
    else
        printf "  %-10s missing\n" "${tool}"
        case "${tool}" in samtools|bedtools|cutadapt) missing=1 ;; esac
    fi
done
for mod in matplotlib numpy; do
    if "${py}" -c "import ${mod}" >/dev/null 2>&1; then
        printf "  %-10s found   (python module, optional - for PDF plots)\n" "${mod}"
    else
        printf "  %-10s missing (python module, optional - for PDF plots)\n" "${mod}"
    fi
done
if ! find_tool STAR >/dev/null && ! find_tool hisat2 >/dev/null; then
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
        echo "STAR ideally needs ${MIN_STAR_RAM_GB} GB RAM for human or mouse genomes. This computer"
        echo "has ${total_gb} GB, so HISAT2 is recommended (run PARALUS with -H). For smaller"
        echo "genomes, STAR also works with less RAM."
        ask "Install HISAT2 instead of STAR?" y && use_hisat=1
    elif [ "${session_gb}" -lt "${MIN_STAR_RAM_GB}" ]; then
        echo
        echo "STAR ideally needs ${MIN_STAR_RAM_GB} GB RAM for human or mouse genomes. This machine has"
        echo "enough, but your current session has ${session_gb} GB, so request more memory for your"
        echo "session or job before running PARALUS with STAR. HISAT2 (-H) needs less"
        echo "memory but reports fewer clusters, so STAR is recommended here."
        ask "Install STAR? Answer n to install HISAT2 instead." y || use_hisat=1
    fi
fi
aligner=$( [ "${use_hisat}" -eq 1 ] && echo hisat2 || echo STAR )
echo "Aligner to install: ${aligner}"

# --------------------------------------------------------------------------
# 4. what will be installed
# --------------------------------------------------------------------------
target=$( [ -n "${prefix}" ] && echo "${prefix}" || echo "conda environment '${ENV_NAME}'" )
update_warning() {
    echo "conda may update some packages already in this environment so that they"
    echo "work with the PARALUS dependencies."
}
echo
line
if [ "${missing}" -eq 0 ]; then
    echo "All required tools are already available, so PARALUS can run without a"
    echo "new environment. A separate environment is still useful to keep the"
    echo "PARALUS tools and versions together."
    if [ "${env_exists}" -eq 1 ]; then
        update_warning
        question="Add all PARALUS dependencies to the existing ${target} anyway?"
    else
        question="Create ${target} with all PARALUS dependencies anyway?"
    fi
    if ! ask "${question}" n; then
        echo
        echo "No environment will be created. You can still set up the paralus command"
        echo "by adding the PARALUS folder to your PATH."
        ask "Set up the paralus command this way?" y || { echo "Stopped."; exit 0; }
        setup_path noask
        exit 0
    fi
elif [ "${env_exists}" -eq 1 ]; then
    echo "The installer will add all PARALUS dependencies to the existing"
    echo "${target}: ${aligner}, samtools, bedtools, cutadapt, python, matplotlib, numpy."
    update_warning
    if ! ask "Continue?" y; then
        # offer a separate environment instead, so nothing in the existing one changes
        if [ -n "${conda_bin}" ] && [ "${ENV_NAME}:${prefix}" != "paralus:" ]; then
            echo
            echo "You can install the PARALUS dependencies into a separate conda environment"
            echo "called 'paralus' instead, so nothing in ${target} changes."
            ask "Use a separate 'paralus' environment?" y || { echo "Stopped."; exit 0; }
            ENV_NAME=paralus; env_given=0; prefix=""
            set_env_path
            target="conda environment '${ENV_NAME}'"
            shortcut_if_complete
            env_exists=0    # the check further down asks again if 'paralus' already exists
        else
            exit 0
        fi
    fi
else
    echo "The installer will create ${target} with all PARALUS"
    echo "dependencies: ${aligner}, samtools, bedtools, cutadapt, python, matplotlib, numpy."
    ask "Continue?" y || exit 0
fi

# STAR already in an existing target environment (e.g. rerun to add new dependencies)
star_in_env=0
[ "${use_hisat}" -eq 0 ] && [ -n "${env_path}" ] && [ -x "${env_path}/bin/STAR" ] && star_in_env=1

# Mac + STAR: reuse an existing compiled STAR, or compile it (needs Homebrew and gcc)
reuse_star=""
if [ "${system}" = mac ] && [ "${use_hisat}" -eq 0 ] && [ "${star_in_env}" -eq 0 ] \
   && existing=$(command -v STAR 2>/dev/null); then
    ver=$(STAR --version 2>/dev/null)
    echo
    echo "Found STAR ${ver} at ${existing}"
    ask "Use this STAR (it will be copied into the new environment)? Answer n to compile a new one." y \
        && reuse_star="${existing}"
fi

if [ "${system}" = mac ] && [ "${use_hisat}" -eq 0 ] && [ "${star_in_env}" -eq 0 ] && [ -z "${reuse_star}" ]; then
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
    if [ "${dry_run}" -eq 1 ] && [ -z "${prefix}" ]; then
        env_flag=(-n "${ENV_NAME}")
        env_path="$HOME/miniconda3/envs/${ENV_NAME}"
        env_desc="the conda environment '${ENV_NAME}'"$'\n'"(${env_path})"
        activate="conda activate ${ENV_NAME}"
    else
        set_env_path
    fi
fi

if [ -d "${env_path}" ]; then
    # already confirmed above when the environment existed from the start
    if [ "${env_exists}" -eq 0 ]; then
        echo
        echo "An environment already exists at ${env_path}."
        update_warning
        ask "Install the PARALUS dependencies into it?" y || { echo "Stopped."; exit 0; }
    fi
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
if [ "${use_hisat}" -eq 0 ] && [ "${star_in_env}" -eq 0 ]; then
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
# 7. the paralus command, check and finish
# --------------------------------------------------------------------------
link_launcher
finish
