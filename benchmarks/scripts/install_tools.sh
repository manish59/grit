#!/usr/bin/env bash
#############################################################################
# Tool Installation and Detection for Multi-Tool Benchmarks
#
# Supports: grit, bedtools, bedops, granges, pyranges, polars-bio
#############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV_DIR="$SCRIPT_DIR/../.venv"

# Colors
if [[ -t 1 ]]; then
  GREEN='\033[0;32m'
  RED='\033[0;31m'
  YELLOW='\033[0;33m'
  CYAN='\033[0;36m'
  NC='\033[0m'
else
  GREEN='' RED='' YELLOW='' CYAN='' NC=''
fi

#############################################################################
# Detection Functions
#############################################################################

detect_grit() {
  if [[ -x "$SCRIPT_DIR/../../target/release/grit" ]]; then
    echo "$SCRIPT_DIR/../../target/release/grit"
  elif command -v grit &>/dev/null; then
    command -v grit
  else
    return 1
  fi
}

detect_bedtools() {
  if command -v bedtools &>/dev/null; then
    command -v bedtools
    return 0
  fi
  return 1
}

detect_bedops() {
  if command -v bedops &>/dev/null; then
    command -v bedops
    return 0
  fi
  return 1
}

detect_granges() {
  if command -v granges &>/dev/null; then
    command -v granges
    return 0
  fi
  return 1
}

detect_python_venv() {
  if [[ -f "$VENV_DIR/bin/python" ]]; then
    echo "$VENV_DIR/bin/python"
    return 0
  fi
  return 1
}

detect_pyranges() {
  local python_bin
  if python_bin=$(detect_python_venv); then
    if "$python_bin" -c "import pyranges" &>/dev/null; then
      return 0
    fi
  fi
  return 1
}

detect_polars_bio() {
  local python_bin
  if python_bin=$(detect_python_venv); then
    if "$python_bin" -c "import polars_bio" &>/dev/null; then
      return 0
    fi
  fi
  return 1
}

#############################################################################
# Version Functions
#############################################################################

get_grit_version() {
  local grit_bin
  if grit_bin=$(detect_grit 2>/dev/null); then
    "$grit_bin" --version 2>/dev/null | head -1 || echo "unknown"
  else
    echo "not installed"
  fi
}

get_bedtools_version() {
  if detect_bedtools &>/dev/null; then
    bedtools --version 2>/dev/null | head -1 || echo "unknown"
  else
    echo "not installed"
  fi
}

get_bedops_version() {
  if detect_bedops &>/dev/null; then
    bedops --version 2>&1 | head -1 || echo "unknown"
  else
    echo "not installed"
  fi
}

get_granges_version() {
  if detect_granges &>/dev/null; then
    granges --version 2>/dev/null | head -1 || echo "unknown"
  else
    echo "not installed"
  fi
}

get_pyranges_version() {
  local python_bin
  if python_bin=$(detect_python_venv 2>/dev/null); then
    "$python_bin" -c "import pyranges; print(f'pyranges {pyranges.__version__}')" 2>/dev/null || echo "not installed"
  else
    echo "not installed"
  fi
}

get_polars_bio_version() {
  local python_bin
  if python_bin=$(detect_python_venv 2>/dev/null); then
    "$python_bin" -c "import polars_bio; print(f'polars-bio {polars_bio.__version__}')" 2>/dev/null || echo "not installed"
  else
    echo "not installed"
  fi
}

#############################################################################
# Installation Functions
#############################################################################

install_bedops() {
  echo -e "${CYAN}Installing bedops...${NC}"

  if [[ "$(uname)" == "Darwin" ]]; then
    if command -v brew &>/dev/null; then
      brew install bedops
    else
      echo -e "${RED}Please install Homebrew first: https://brew.sh${NC}"
      return 1
    fi
  else
    # Linux - try apt, yum, or conda
    if command -v apt-get &>/dev/null; then
      sudo apt-get update && sudo apt-get install -y bedops
    elif command -v yum &>/dev/null; then
      sudo yum install -y bedops
    elif command -v conda &>/dev/null; then
      conda install -c bioconda bedops
    else
      echo -e "${YELLOW}Please install bedops manually: https://bedops.readthedocs.io/en/latest/content/installation.html${NC}"
      return 1
    fi
  fi

  echo -e "${GREEN}bedops installed successfully${NC}"
}

install_granges() {
  echo -e "${CYAN}Installing granges...${NC}"

  if command -v cargo &>/dev/null; then
    cargo install granges
    echo -e "${GREEN}granges installed successfully${NC}"
  else
    echo -e "${RED}Please install Rust first: https://rustup.rs${NC}"
    return 1
  fi
}

setup_python_venv() {
  echo -e "${CYAN}Setting up Python virtual environment...${NC}"

  if [[ ! -d "$VENV_DIR" ]]; then
    python3 -m venv "$VENV_DIR"
  fi

  "$VENV_DIR/bin/pip" install --upgrade pip
  echo -e "${GREEN}Python venv ready at $VENV_DIR${NC}"
}

install_pyranges() {
  echo -e "${CYAN}Installing pyranges...${NC}"

  setup_python_venv
  "$VENV_DIR/bin/pip" install pyranges pandas

  echo -e "${GREEN}pyranges installed successfully${NC}"
}

install_polars_bio() {
  echo -e "${CYAN}Installing polars-bio...${NC}"

  setup_python_venv
  "$VENV_DIR/bin/pip" install polars-bio polars

  echo -e "${GREEN}polars-bio installed successfully${NC}"
}

install_all_python() {
  setup_python_venv
  "$VENV_DIR/bin/pip" install pyranges polars-bio polars pandas
  echo -e "${GREEN}All Python packages installed${NC}"
}

#############################################################################
# Status Functions
#############################################################################

print_status() {
  echo -e "${CYAN}=== Benchmark Tool Status ===${NC}"
  echo ""

  # CLI Tools
  echo -e "${CYAN}CLI Tools:${NC}"

  printf "  %-12s " "grit:"
  if detect_grit &>/dev/null; then
    echo -e "${GREEN}$(get_grit_version)${NC}"
  else
    echo -e "${RED}not installed${NC} (build with: cargo build --release)"
  fi

  printf "  %-12s " "bedtools:"
  if detect_bedtools &>/dev/null; then
    echo -e "${GREEN}$(get_bedtools_version)${NC}"
  else
    echo -e "${RED}not installed${NC}"
  fi

  printf "  %-12s " "bedops:"
  if detect_bedops &>/dev/null; then
    echo -e "${GREEN}$(get_bedops_version)${NC}"
  else
    echo -e "${YELLOW}not installed${NC} (run: ./install_tools.sh install bedops)"
  fi

  printf "  %-12s " "granges:"
  if detect_granges &>/dev/null; then
    echo -e "${GREEN}$(get_granges_version)${NC}"
  else
    echo -e "${YELLOW}not installed${NC} (run: ./install_tools.sh install granges)"
  fi

  echo ""
  echo -e "${CYAN}Python Libraries:${NC}"

  printf "  %-12s " "pyranges:"
  echo -e "$(get_pyranges_version)"

  printf "  %-12s " "polars-bio:"
  echo -e "$(get_polars_bio_version)"

  echo ""
}

print_available_tools() {
  local available=""

  detect_grit &>/dev/null && available="$available grit"
  detect_bedtools &>/dev/null && available="$available bedtools"
  detect_bedops &>/dev/null && available="$available bedops"
  detect_granges &>/dev/null && available="$available granges"
  detect_pyranges && available="$available pyranges"
  detect_polars_bio && available="$available polars-bio"

  echo "$available" | xargs
}

#############################################################################
# Main
#############################################################################

usage() {
  cat <<EOF
Benchmark Tool Installation & Detection

Usage:
  $0 status              Show installation status of all tools
  $0 available           List available (installed) tools
  $0 install <tool>      Install a specific tool
  $0 install-python      Install all Python packages (pyranges, polars-bio)
  $0 detect <tool>       Check if tool is installed (exit 0/1)
  $0 version <tool>      Get version of installed tool

Tools: bedops, granges, pyranges, polars-bio

Examples:
  $0 status
  $0 install bedops
  $0 install-python
  $0 detect pyranges && echo "pyranges available"
EOF
}

main() {
  local cmd="${1:-status}"

  case "$cmd" in
    status)
      print_status
      ;;
    available)
      print_available_tools
      ;;
    detect)
      local tool="${2:-}"
      case "$tool" in
        grit) detect_grit >/dev/null ;;
        bedtools) detect_bedtools >/dev/null ;;
        bedops) detect_bedops >/dev/null ;;
        granges) detect_granges >/dev/null ;;
        pyranges) detect_pyranges ;;
        polars-bio|polars_bio) detect_polars_bio ;;
        *) echo "Unknown tool: $tool"; exit 1 ;;
      esac
      ;;
    version)
      local tool="${2:-}"
      case "$tool" in
        grit) get_grit_version ;;
        bedtools) get_bedtools_version ;;
        bedops) get_bedops_version ;;
        granges) get_granges_version ;;
        pyranges) get_pyranges_version ;;
        polars-bio|polars_bio) get_polars_bio_version ;;
        *) echo "Unknown tool: $tool"; exit 1 ;;
      esac
      ;;
    install)
      local tool="${2:-}"
      case "$tool" in
        bedops) install_bedops ;;
        granges) install_granges ;;
        pyranges) install_pyranges ;;
        polars-bio|polars_bio) install_polars_bio ;;
        *) echo "Unknown tool: $tool"; usage; exit 1 ;;
      esac
      ;;
    install-python)
      install_all_python
      ;;
    help|--help|-h)
      usage
      ;;
    *)
      usage
      exit 1
      ;;
  esac
}

# Export functions for sourcing
export -f detect_grit detect_bedtools detect_bedops detect_granges
export -f detect_python_venv detect_pyranges detect_polars_bio

# Only run main if executed directly (not sourced)
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  main "$@"
fi
