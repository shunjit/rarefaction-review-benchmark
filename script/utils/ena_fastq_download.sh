#!/bin/bash
#===============================================================================
# ENA/NCBI/DDBJ Universal FASTQ Downloader
#===============================================================================
# 概要:
#   INSDC (ENA/NCBI SRA/DDBJ) からFASTQファイルを並列ダウンロードし、
#   MD5チェックサムで検証するスクリプト
#
# 対応アクセション:
#   - BioProject: PRJNA*, PRJEB*, PRJD* (プロジェクト全体)
#   - Study:      SRP*, ERP*, DRP*      (スタディ単位)
#   - Sample:     SRS*, ERS*, DRS*      (サンプル単位)
#   - Experiment: SRX*, ERX*, DRX*      (実験単位)
#   - Run:        SRR*, ERR*, DRR*      (個別Run)
#
# 使用方法:
#   ./ena_fastq_download.sh <ACCESSION> [OPTIONS]
#
# 必須引数:
#   ACCESSION    : アクセション番号（例: PRJNA255789, SRR12345678）
#
# オプション引数:
#   -o, --output DIR          : 出力ベースディレクトリ（デフォルト: カレント）
#   -p, --parallel N          : 並列ダウンロード数（デフォルト: 2）
#   -f, --filter-strategy STR : library_strategyでフィルター（例: AMPLICON, WGS）
#   -P, --filter-platform STR : instrument_platformでフィルター（例: ILLUMINA, LS454）
#   -m, --md5                 : MD5検証を実施（デフォルト）
#   -M, --no-md5              : MD5検証をスキップ
#   -s, --skip-existing       : 既存ファイルをスキップ（デフォルト）
#   -S, --no-skip             : 既存ファイルを再ダウンロード
#   -d, --dry-run             : ダウンロードせずメタデータのみ取得
#   -h, --help                : ヘルプを表示
#
# 使用例:
#   # NCBIプロジェクトをダウンロード
#   ./ena_fastq_download.sh PRJNA255789
#
#   # EBIプロジェクトを指定ディレクトリへダウンロード
#   ./ena_fastq_download.sh PRJEB10725 -o /Volumes/PS3000/seq_download
#
#   # 並列数4、MD5スキップで高速ダウンロード
#   ./ena_fastq_download.sh PRJNA255789 -p 4 --no-md5
#
#   # 特定のRunのみダウンロード
#   ./ena_fastq_download.sh SRR1234567
#
#   # ドライラン（メタデータ確認のみ）
#   ./ena_fastq_download.sh PRJNA255789 --dry-run
#
#   # 16S rRNAアンプリコンデータのみダウンロード
#   ./ena_fastq_download.sh PRJEB10725 --filter-strategy AMPLICON
#
#   # メタゲノムWGSデータのみダウンロード
#   ./ena_fastq_download.sh PRJEB10725 -f WGS
#
#   # Illuminaデータのみダウンロード
#   ./ena_fastq_download.sh PRJEB10725 --filter-platform ILLUMINA
#
#   # フィルター組み合わせ: Illumina 16S rRNAのみ
#   ./ena_fastq_download.sh PRJEB10725 -f AMPLICON -P ILLUMINA
#
# 動作環境:
#   macOS (Sequoia 15.x, Apple Silicon) / Linux
#   必須: curl, wget, awk, md5 (macOS) or md5sum (Linux)
#
# データソース:
#   ENA Portal API (https://www.ebi.ac.uk/ena/portal/api/)
#   INSDCデータ共有により、NCBI/DDBJのデータもENA経由で取得可能
#
# 作成日: 2025-01
# バージョン: 2.2.0
# 変更履歴:
#   2.2.0 - --filter-platform オプション追加（instrument_platformでフィルタリング）
#         - ドライラン時にプラットフォーム別・クロス集計を表示
#         - フィルターの組み合わせに対応（-f と -P の同時使用）
#   2.1.0 - --filter-strategy オプション追加（library_strategyでフィルタリング）
#         - ドライラン時にストラテジー別統計を表示
#   2.0.0 - 初版
#===============================================================================

set -euo pipefail

#-------------------------------------------------------------------------------
# バージョン情報
#-------------------------------------------------------------------------------
readonly SCRIPT_NAME="ena_fastq_download.sh"
readonly SCRIPT_VERSION="2.2.0"

#-------------------------------------------------------------------------------
# デフォルト設定
#-------------------------------------------------------------------------------
DEFAULT_OUTPUT_DIR="."
DEFAULT_PARALLEL=2
DEFAULT_VERIFY_MD5=1
DEFAULT_SKIP_EXISTING=1
DEFAULT_DRY_RUN=0
DEFAULT_FILTER_STRATEGY=""  # 空=全て、WGS/AMPLICON等で絞り込み
DEFAULT_FILTER_PLATFORM=""  # 空=全て、ILLUMINA/LS454等で絞り込み

# wget設定
WGET_RETRIES=3
WGET_WAIT=5
WGET_TIMEOUT=60
WGET_READ_TIMEOUT=300

# ディスク容量マージン
DISK_MARGIN=1.2

#-------------------------------------------------------------------------------
# グローバル変数（引数解析で設定）
#-------------------------------------------------------------------------------
ACCESSION=""
OUTPUT_DIR="${DEFAULT_OUTPUT_DIR}"
PARALLEL_DL="${DEFAULT_PARALLEL}"
VERIFY_MD5="${DEFAULT_VERIFY_MD5}"
SKIP_EXISTING="${DEFAULT_SKIP_EXISTING}"
DRY_RUN="${DEFAULT_DRY_RUN}"
FILTER_STRATEGY="${DEFAULT_FILTER_STRATEGY}"
FILTER_PLATFORM="${DEFAULT_FILTER_PLATFORM}"

#-------------------------------------------------------------------------------
# カラー出力設定
#-------------------------------------------------------------------------------
setup_colors() {
    if [[ -t 1 ]]; then
        COLOR_RED='\033[0;31m'
        COLOR_GREEN='\033[0;32m'
        COLOR_YELLOW='\033[0;33m'
        COLOR_BLUE='\033[0;34m'
        COLOR_CYAN='\033[0;36m'
        COLOR_BOLD='\033[1m'
        COLOR_RESET='\033[0m'
    else
        COLOR_RED=''
        COLOR_GREEN=''
        COLOR_YELLOW=''
        COLOR_BLUE=''
        COLOR_CYAN=''
        COLOR_BOLD=''
        COLOR_RESET=''
    fi
}

#-------------------------------------------------------------------------------
# ログ関数
#-------------------------------------------------------------------------------
log_info() {
    echo -e "${COLOR_BLUE}[$(date '+%H:%M:%S')]${COLOR_RESET} $*"
}

log_success() {
    echo -e "${COLOR_GREEN}[$(date '+%H:%M:%S')]${COLOR_RESET} ✓ $*"
}

log_warn() {
    echo -e "${COLOR_YELLOW}[$(date '+%H:%M:%S')]${COLOR_RESET} ⚠ $*" >&2
}

log_error() {
    echo -e "${COLOR_RED}[$(date '+%H:%M:%S')]${COLOR_RESET} ✗ $*" >&2
}

log_header() {
    echo ""
    echo -e "${COLOR_BOLD}${COLOR_CYAN}$*${COLOR_RESET}"
    echo -e "${COLOR_CYAN}$(printf '%.0s─' {1..50})${COLOR_RESET}"
}

#-------------------------------------------------------------------------------
# ヘルプ表示
#-------------------------------------------------------------------------------
show_help() {
    cat << EOF
${SCRIPT_NAME} v${SCRIPT_VERSION}
ENA/NCBI/DDBJ Universal FASTQ Downloader

USAGE:
    ${SCRIPT_NAME} <ACCESSION> [OPTIONS]

ARGUMENTS:
    ACCESSION           Accession number (PRJNA*, PRJEB*, SRR*, etc.)

OPTIONS:
    -o, --output DIR          Output base directory (default: current directory)
    -p, --parallel N          Number of parallel downloads (default: 2)
    -f, --filter-strategy STR Filter by library_strategy (e.g., AMPLICON, WGS)
    -P, --filter-platform STR Filter by instrument_platform (e.g., ILLUMINA, LS454)
    -m, --md5                 Verify MD5 checksums (default)
    -M, --no-md5              Skip MD5 verification (faster)
    -s, --skip-existing       Skip existing files with matching size (default)
    -S, --no-skip             Re-download existing files
    -d, --dry-run             Fetch metadata only, don't download
    -h, --help                Show this help message

SUPPORTED ACCESSION TYPES:
    BioProject    PRJNA*, PRJEB*, PRJD*   (entire project)
    Study         SRP*, ERP*, DRP*        (study level)
    Sample        SRS*, ERS*, DRS*        (sample level)
    Experiment    SRX*, ERX*, DRX*        (experiment level)
    Run           SRR*, ERR*, DRR*        (single run)

LIBRARY STRATEGIES (for --filter-strategy):
    AMPLICON      Amplicon sequencing (e.g., 16S rRNA)
    WGS           Whole Genome Shotgun (metagenomic)
    WXS           Whole Exome Sequencing
    RNA-Seq       RNA sequencing
    (Use --dry-run first to see available strategies in a project)

INSTRUMENT PLATFORMS (for --filter-platform):
    ILLUMINA      Illumina platforms (HiSeq, MiSeq, NovaSeq, etc.)
    LS454         Roche 454 pyrosequencing
    OXFORD_NANOPORE  Oxford Nanopore long-read
    PACBIO_SMRT   PacBio long-read
    ION_TORRENT   Ion Torrent platforms
    (Use --dry-run first to see available platforms in a project)

EXAMPLES:
    # Download NCBI project
    ${SCRIPT_NAME} PRJNA255789

    # Download to specific directory with 4 parallel connections
    ${SCRIPT_NAME} PRJEB10725 -o /Volumes/SSD/data -p 4

    # Fast download (skip MD5 verification)
    ${SCRIPT_NAME} PRJNA255789 --no-md5

    # Download specific run only
    ${SCRIPT_NAME} SRR1234567

    # Check project size before downloading (recommended first step)
    ${SCRIPT_NAME} PRJNA255789 --dry-run

    # Download only 16S rRNA amplicon data (filter by strategy)
    ${SCRIPT_NAME} PRJEB10725 --filter-strategy AMPLICON

    # Download only metagenomic WGS data
    ${SCRIPT_NAME} PRJEB10725 -f WGS

    # Download only Illumina data (filter by platform)
    ${SCRIPT_NAME} PRJEB10725 --filter-platform ILLUMINA

    # Combine filters: Illumina 16S rRNA amplicon only
    ${SCRIPT_NAME} PRJEB10725 -f AMPLICON -P ILLUMINA

EOF
    exit 0
}

#-------------------------------------------------------------------------------
# 引数解析
#-------------------------------------------------------------------------------
parse_args() {
    if [[ $# -eq 0 ]]; then
        log_error "Accession number is required"
        echo "Usage: ${SCRIPT_NAME} <ACCESSION> [OPTIONS]"
        echo "Run '${SCRIPT_NAME} --help' for more information"
        exit 1
    fi

    while [[ $# -gt 0 ]]; do
        case "$1" in
            -h|--help)
                show_help
                ;;
            -o|--output)
                OUTPUT_DIR="$2"
                shift 2
                ;;
            -p|--parallel)
                PARALLEL_DL="$2"
                shift 2
                ;;
            -m|--md5)
                VERIFY_MD5=1
                shift
                ;;
            -M|--no-md5)
                VERIFY_MD5=0
                shift
                ;;
            -s|--skip-existing)
                SKIP_EXISTING=1
                shift
                ;;
            -S|--no-skip)
                SKIP_EXISTING=0
                shift
                ;;
            -d|--dry-run)
                DRY_RUN=1
                shift
                ;;
            -f|--filter-strategy)
                FILTER_STRATEGY=$(echo "$2" | tr '[:lower:]' '[:upper:]')
                shift 2
                ;;
            -P|--filter-platform)
                FILTER_PLATFORM=$(echo "$2" | tr '[:lower:]' '[:upper:]')
                shift 2
                ;;
            -*)
                log_error "Unknown option: $1"
                exit 1
                ;;
            *)
                if [[ -z "${ACCESSION}" ]]; then
                    ACCESSION="$1"
                else
                    log_error "Multiple accessions not supported. Use one accession per run."
                    exit 1
                fi
                shift
                ;;
        esac
    done

    # アクセション必須チェック
    if [[ -z "${ACCESSION}" ]]; then
        log_error "Accession number is required"
        exit 1
    fi
}

#-------------------------------------------------------------------------------
# アクセション形式の検証と分類
#-------------------------------------------------------------------------------
validate_accession() {
    local acc="$1"
    local acc_upper
    acc_upper=$(echo "${acc}" | tr '[:lower:]' '[:upper:]')
    
    # 対応するアクセション形式のパターン
    # BioProject: PRJNA, PRJEB, PRJD, PRJDA, PRJDB
    # Study: SRP, ERP, DRP
    # Sample: SRS, ERS, DRS, SAMN, SAME, SAMD
    # Experiment: SRX, ERX, DRX
    # Run: SRR, ERR, DRR
    
    if [[ "${acc_upper}" =~ ^PRJ[NEDB][A-Z]?[0-9]+$ ]]; then
        echo "bioproject"
    elif [[ "${acc_upper}" =~ ^[SED]RP[0-9]+$ ]]; then
        echo "study"
    elif [[ "${acc_upper}" =~ ^[SED]RS[0-9]+$ ]] || [[ "${acc_upper}" =~ ^SAM[NED][A-Z]?[0-9]+$ ]]; then
        echo "sample"
    elif [[ "${acc_upper}" =~ ^[SED]RX[0-9]+$ ]]; then
        echo "experiment"
    elif [[ "${acc_upper}" =~ ^[SED]RR[0-9]+$ ]]; then
        echo "run"
    else
        echo "unknown"
    fi
}

#-------------------------------------------------------------------------------
# ユーティリティ関数
#-------------------------------------------------------------------------------

# バイト数を人間可読形式に変換
human_bytes() {
    local bytes="$1"
    awk -v b="${bytes}" 'BEGIN {
        units[0]="B"; units[1]="KB"; units[2]="MB"; units[3]="GB"; units[4]="TB"
        for(i=0; b>=1024 && i<4; i++) b/=1024
        printf "%.2f %s\n", b, units[i]
    }'
}

# 秒数を人間可読形式に変換
human_time() {
    local seconds="$1"
    if (( seconds < 60 )); then
        echo "${seconds}s"
    elif (( seconds < 3600 )); then
        printf "%dm %ds" $((seconds / 60)) $((seconds % 60))
    else
        printf "%dh %dm %ds" $((seconds / 3600)) $((seconds % 3600 / 60)) $((seconds % 60))
    fi
}

# OS検出
detect_os() {
    case "$(uname -s)" in
        Darwin) echo "macos" ;;
        Linux)  echo "linux" ;;
        *)      echo "unknown" ;;
    esac
}

# MD5計算（OS互換）
calculate_md5() {
    local filepath="$1"
    local os
    os=$(detect_os)
    
    if [[ "${os}" == "macos" ]]; then
        md5 -q "${filepath}"
    else
        md5sum "${filepath}" | awk '{print $1}'
    fi
}

# ファイルサイズ取得（OS互換）
get_file_size() {
    local filepath="$1"
    local os
    os=$(detect_os)
    
    if [[ "${os}" == "macos" ]]; then
        stat -f%z "${filepath}" 2>/dev/null || echo "0"
    else
        stat --printf="%s" "${filepath}" 2>/dev/null || echo "0"
    fi
}

# 必須コマンドの確認
check_dependencies() {
    local missing=()
    local os
    os=$(detect_os)
    
    for cmd in curl wget awk; do
        if ! command -v "${cmd}" &>/dev/null; then
            missing+=("${cmd}")
        fi
    done
    
    # MD5コマンド（OS依存）
    if [[ "${os}" == "macos" ]]; then
        if ! command -v md5 &>/dev/null; then
            missing+=("md5")
        fi
    else
        if ! command -v md5sum &>/dev/null; then
            missing+=("md5sum")
        fi
    fi
    
    if [[ ${#missing[@]} -gt 0 ]]; then
        log_error "Missing required commands: ${missing[*]}"
        exit 1
    fi
}

# ディスク空き容量チェック
check_disk_space() {
    local required_bytes="$1"
    local target_dir="$2"
    
    local available_bytes
    available_bytes=$(df -P "${target_dir}" | awk 'NR==2 {print $4 * 1024}')
    
    local required_with_margin
    required_with_margin=$(awk -v r="${required_bytes}" -v m="${DISK_MARGIN}" \
        'BEGIN {printf "%.0f", r * m}')
    
    log_info "Required (×${DISK_MARGIN}): $(human_bytes "${required_with_margin}")"
    log_info "Available: $(human_bytes "${available_bytes}")"
    
    if (( available_bytes < required_with_margin )); then
        log_error "Insufficient disk space!"
        return 1
    fi
    
    log_success "Disk space OK"
    return 0
}

#-------------------------------------------------------------------------------
# ENA API関数
#-------------------------------------------------------------------------------

# ENA APIからメタデータ取得
fetch_ena_metadata() {
    local accession="$1"
    local output_file="$2"
    
    local api_url="https://www.ebi.ac.uk/ena/portal/api/filereport"
    api_url+="?accession=${accession}"
    api_url+="&result=read_run"
    api_url+="&fields=run_accession,sample_accession,experiment_accession,study_accession"
    api_url+=",library_layout,library_source,library_strategy,library_selection"
    api_url+=",instrument_platform,instrument_model"
    api_url+=",fastq_ftp,fastq_md5,fastq_bytes"
    api_url+=",read_count,base_count"
    api_url+="&format=tsv"
    api_url+="&download=true"
    api_url+="&limit=0"
    
    if ! curl -fsSL "${api_url}" > "${output_file}" 2>/dev/null; then
        return 1
    fi
    
    # データが取得できたか確認（ヘッダー行のみでないか）
    local line_count
    line_count=$(wc -l < "${output_file}" | tr -d ' ')
    
    if [[ ${line_count} -le 1 ]]; then
        return 1
    fi
    
    return 0
}

# メタデータからURL/MD5/サイズを展開（フィルター適用）
parse_metadata() {
    local meta_file="$1"
    local output_file="$2"
    local filter_strategy="${3:-}"  # オプション: フィルターストラテジー
    local filter_platform="${4:-}"  # オプション: フィルタープラットフォーム
    
    awk -F'\t' -v filter_str="${filter_strategy}" -v filter_plat="${filter_platform}" '
    NR==1 {
        for(i=1; i<=NF; i++) { h[$i]=i }
        next
    }
    {
        # library_strategyフィルター適用
        if(filter_str != "" && $(h["library_strategy"]) != filter_str) next
        
        # instrument_platformフィルター適用
        if(filter_plat != "" && $(h["instrument_platform"]) != filter_plat) next
        
        ftp = $(h["fastq_ftp"])
        md5 = $(h["fastq_md5"])
        byt = $(h["fastq_bytes"])
        run = $(h["run_accession"])
        
        if(ftp == "") next
        
        n = split(ftp, urls, ";")
        split(md5, md5s, ";")
        split(byt, sizes, ";")
        
        for(i=1; i<=n; i++) {
            if(urls[i] != "") {
                # HTTPSスキームを付与
                print "https://" urls[i] "\t" md5s[i] "\t" sizes[i] "\t" run
            }
        }
    }
    ' "${meta_file}" > "${output_file}"
}

# library_strategy別の統計を表示
show_strategy_breakdown() {
    local meta_file="$1"
    
    echo ""
    log_info "Breakdown by library_strategy:"
    
    awk -F'\t' '
    NR==1 {
        for(i=1; i<=NF; i++) { h[$i]=i }
        next
    }
    {
        strategy = $(h["library_strategy"])
        n = split($(h["fastq_bytes"]), bytes, ";")
        total = 0
        for(i=1; i<=n; i++) { total += bytes[i]; files++ }
        
        count[strategy]++
        size[strategy] += total
        filecount[strategy] += n
    }
    END {
        # サイズでソート（降順）
        for(s in size) {
            printf "%012.0f\t%s\t%d\t%d\n", size[s], s, count[s], filecount[s]
        }
    }
    ' "${meta_file}" | sort -rn | while IFS=$'\t' read -r bytes strategy runs files; do
        local human_size
        human_size=$(human_bytes "${bytes}")
        printf "  ${COLOR_CYAN}%-12s${COLOR_RESET}: %3d runs, %3d files, %s\n" \
            "${strategy}" "${runs}" "${files}" "${human_size}"
    done
}

# instrument_platform別の統計を表示
show_platform_breakdown() {
    local meta_file="$1"
    
    echo ""
    log_info "Breakdown by instrument_platform:"
    
    awk -F'\t' '
    NR==1 {
        for(i=1; i<=NF; i++) { h[$i]=i }
        next
    }
    {
        platform = $(h["instrument_platform"])
        n = split($(h["fastq_bytes"]), bytes, ";")
        total = 0
        for(i=1; i<=n; i++) { total += bytes[i]; files++ }
        
        count[platform]++
        size[platform] += total
        filecount[platform] += n
    }
    END {
        # サイズでソート（降順）
        for(p in size) {
            printf "%012.0f\t%s\t%d\t%d\n", size[p], p, count[p], filecount[p]
        }
    }
    ' "${meta_file}" | sort -rn | while IFS=$'\t' read -r bytes platform runs files; do
        local human_size
        human_size=$(human_bytes "${bytes}")
        printf "  ${COLOR_CYAN}%-16s${COLOR_RESET}: %3d runs, %3d files, %s\n" \
            "${platform}" "${runs}" "${files}" "${human_size}"
    done
}

# Strategy × Platform のクロス集計を表示
show_cross_breakdown() {
    local meta_file="$1"
    
    echo ""
    log_info "Breakdown by library_strategy × instrument_platform:"
    
    awk -F'\t' '
    NR==1 {
        for(i=1; i<=NF; i++) { h[$i]=i }
        next
    }
    {
        strategy = $(h["library_strategy"])
        platform = $(h["instrument_platform"])
        key = strategy "\t" platform
        
        n = split($(h["fastq_bytes"]), bytes, ";")
        total = 0
        for(i=1; i<=n; i++) { total += bytes[i]; files++ }
        
        count[key]++
        size[key] += total
        filecount[key] += n
    }
    END {
        for(k in size) {
            printf "%012.0f\t%s\t%d\t%d\n", size[k], k, count[k], filecount[k]
        }
    }
    ' "${meta_file}" | sort -rn | while IFS=$'\t' read -r bytes strategy platform runs files; do
        local human_size
        human_size=$(human_bytes "${bytes}")
        printf "  ${COLOR_CYAN}%-12s${COLOR_RESET} × ${COLOR_CYAN}%-16s${COLOR_RESET}: %3d runs, %3d files, %s\n" \
            "${strategy}" "${platform}" "${runs}" "${files}" "${human_size}"
    done
}

# メタデータサマリー表示
show_metadata_summary() {
    local meta_file="$1"
    
    # Run数カウント
    local run_count
    run_count=$(($(wc -l < "${meta_file}") - 1))
    
    # ライブラリレイアウト集計
    local layout_summary
    layout_summary=$(awk -F'\t' 'NR>1 {print $5}' "${meta_file}" | sort | uniq -c | \
        awk '{printf "%s: %d, ", $2, $1}' | sed 's/, $//')
    
    # プラットフォーム集計
    local platform_summary
    platform_summary=$(awk -F'\t' 'NR>1 && NF>=10 {print $10}' "${meta_file}" | sort | uniq -c | \
        awk '{printf "%s: %d, ", $2, $1}' | sed 's/, $//')
    
    # リードカウント合計
    local total_reads
    total_reads=$(awk -F'\t' 'NR>1 {s+=$NF} END{printf "%.0f", s}' "${meta_file}" 2>/dev/null || echo "N/A")
    
    log_info "Runs: ${run_count}"
    log_info "Layout: ${layout_summary:-N/A}"
    log_info "Platform: ${platform_summary:-N/A}"
    if [[ "${total_reads}" != "N/A" && "${total_reads}" != "0" ]]; then
        log_info "Total reads: ${total_reads}"
    fi
}

#-------------------------------------------------------------------------------
# メイン処理
#-------------------------------------------------------------------------------

main() {
    local script_start
    script_start=$(date +%s)
    
    # カラー設定
    setup_colors
    
    # 引数解析
    parse_args "$@"
    
    # アクセション検証
    local acc_type
    acc_type=$(validate_accession "${ACCESSION}")
    
    if [[ "${acc_type}" == "unknown" ]]; then
        log_error "Unknown or unsupported accession format: ${ACCESSION}"
        log_error "Supported formats: PRJNA*, PRJEB*, PRJD*, SRP*, SRR*, etc."
        exit 1
    fi
    
    # 作業ディレクトリ設定
    local work_dir="${OUTPUT_DIR}/${ACCESSION}"
    
    # ヘッダー表示
    echo ""
    echo -e "${COLOR_BOLD}╔══════════════════════════════════════════════════╗${COLOR_RESET}"
    echo -e "${COLOR_BOLD}║  ENA/NCBI/DDBJ FASTQ Downloader v${SCRIPT_VERSION}          ║${COLOR_RESET}"
    echo -e "${COLOR_BOLD}╚══════════════════════════════════════════════════╝${COLOR_RESET}"
    echo ""
    log_info "Accession:     ${ACCESSION} (${acc_type})"
    log_info "Output:        ${work_dir}"
    log_info "Parallel:      ${PARALLEL_DL}"
    log_info "Verify MD5:    $([ "${VERIFY_MD5}" -eq 1 ] && echo 'Yes' || echo 'No')"
    log_info "Skip existing: $([ "${SKIP_EXISTING}" -eq 1 ] && echo 'Yes' || echo 'No')"
    if [[ -n "${FILTER_STRATEGY}" ]]; then
        log_info "Filter:        library_strategy=${FILTER_STRATEGY}"
    fi
    if [[ -n "${FILTER_PLATFORM}" ]]; then
        log_info "Filter:        instrument_platform=${FILTER_PLATFORM}"
    fi
    if [[ "${DRY_RUN}" -eq 1 ]]; then
        log_warn "DRY RUN MODE - No files will be downloaded"
    fi
    
    #---------------------------------------------------------------------------
    # 前提条件チェック
    #---------------------------------------------------------------------------
    log_header "Checking Prerequisites"
    check_dependencies
    log_success "All dependencies available"
    
    #---------------------------------------------------------------------------
    # 作業ディレクトリ準備
    #---------------------------------------------------------------------------
    log_header "Preparing Directories"
    mkdir -p "${work_dir}"/{00_meta,01_fastq,99_logs}
    
    cd "${work_dir}" || {
        log_error "Cannot change to: ${work_dir}"
        exit 1
    }
    log_info "Working directory: $(pwd)"
    
    #---------------------------------------------------------------------------
    # Step 1: メタデータ取得
    #---------------------------------------------------------------------------
    log_header "Step 1: Fetching Metadata from ENA"
    
    local meta_file="00_meta/${ACCESSION}.metadata.tsv"
    
    if ! fetch_ena_metadata "${ACCESSION}" "${meta_file}"; then
        log_error "Failed to fetch metadata for: ${ACCESSION}"
        log_error ""
        log_error "Possible causes:"
        log_error "  1. Invalid accession number"
        log_error "  2. Data not yet synced to ENA (NCBI data may take 24-48h)"
        log_error "  3. No FASTQ files available (raw data only)"
        log_error "  4. Network/API issues"
        log_error ""
        log_error "Try verifying at: https://www.ebi.ac.uk/ena/browser/view/${ACCESSION}"
        exit 1
    fi
    
    log_success "Metadata retrieved"
    show_metadata_summary "${meta_file}"
    
    #---------------------------------------------------------------------------
    # Step 2: ダウンロード情報展開
    #---------------------------------------------------------------------------
    log_header "Step 2: Parsing Download Information"
    
    # フィルターなしの統計を先に表示（ドライラン時）
    if [[ "${DRY_RUN}" -eq 1 ]]; then
        show_strategy_breakdown "${meta_file}"
        show_platform_breakdown "${meta_file}"
        show_cross_breakdown "${meta_file}"
    fi
    
    # フィルター適用してURL/MD5/サイズを展開
    parse_metadata "${meta_file}" "00_meta/${ACCESSION}.urls.md5.bytes.tsv" "${FILTER_STRATEGY}" "${FILTER_PLATFORM}"
    cut -f1 "00_meta/${ACCESSION}.urls.md5.bytes.tsv" > "00_meta/${ACCESSION}.urls.txt"
    
    local file_count
    file_count=$(wc -l < "00_meta/${ACCESSION}.urls.txt" | tr -d ' ')
    
    local total_bytes
    total_bytes=$(awk -F'\t' '{s+=$3} END{print s}' "00_meta/${ACCESSION}.urls.md5.bytes.tsv")
    echo "${total_bytes}" > "99_logs/total_bytes.txt"
    
    log_info "FASTQ files: ${file_count}"
    log_info "Total size:  $(human_bytes "${total_bytes}")"
    
    # メタデータをログに保存
    cp "${meta_file}" "99_logs/"
    
    #---------------------------------------------------------------------------
    # ドライランの場合はここで終了
    #---------------------------------------------------------------------------
    if [[ "${DRY_RUN}" -eq 1 ]]; then
        log_header "Dry Run Summary"
        log_info "Accession:   ${ACCESSION}"
        if [[ -n "${FILTER_STRATEGY}" ]]; then
            log_info "Filter:      library_strategy=${FILTER_STRATEGY}"
        fi
        if [[ -n "${FILTER_PLATFORM}" ]]; then
            log_info "Filter:      instrument_platform=${FILTER_PLATFORM}"
        fi
        log_info "Files:       ${file_count}"
        log_info "Total size:  $(human_bytes "${total_bytes}")"
        log_info ""
        log_info "Metadata saved to: ${work_dir}/00_meta/"
        log_info "To download, run without --dry-run option"
        
        # フィルターが未設定で複数ストラテジー/プラットフォームがある場合、ヒントを表示
        if [[ -z "${FILTER_STRATEGY}" && -z "${FILTER_PLATFORM}" ]]; then
            local strategy_count platform_count
            strategy_count=$(awk -F'\t' 'NR>1 {print $7}' "${meta_file}" | sort -u | wc -l | tr -d ' ')
            platform_count=$(awk -F'\t' 'NR>1 {print $9}' "${meta_file}" | sort -u | wc -l | tr -d ' ')
            
            if [[ ${strategy_count} -gt 1 || ${platform_count} -gt 1 ]]; then
                log_info ""
                log_info "💡 TIP: This project contains multiple data types."
                if [[ ${strategy_count} -gt 1 ]]; then
                    log_info "   Filter by library_strategy:"
                    log_info "     -f AMPLICON  (for 16S rRNA amplicon)"
                    log_info "     -f WGS       (for metagenomic shotgun)"
                fi
                if [[ ${platform_count} -gt 1 ]]; then
                    log_info "   Filter by instrument_platform:"
                    log_info "     -P ILLUMINA  (for Illumina data only)"
                    log_info "     -P LS454     (for 454 pyrosequencing)"
                fi
                log_info "   Combine filters for precise selection:"
                log_info "     -f AMPLICON -P ILLUMINA"
            fi
        fi
        exit 0
    fi
    
    #---------------------------------------------------------------------------
    # Step 3: ディスク容量チェック
    #---------------------------------------------------------------------------
    log_header "Step 3: Checking Disk Space"
    
    if ! check_disk_space "${total_bytes}" "${work_dir}"; then
        exit 1
    fi
    
    #---------------------------------------------------------------------------
    # Step 4: ダウンロード
    #---------------------------------------------------------------------------
    log_header "Step 4: Downloading Files"
    
    local dl_start
    dl_start=$(date +%s)
    
    # スキップ対象を除外したダウンロードリスト作成
    local skipped=0
    : > "00_meta/${ACCESSION}.urls_to_download.txt"
    
    while IFS=$'\t' read -r url md5 expected_size run_acc; do
        local filename
        filename="$(basename "${url}")"
        local filepath="01_fastq/${filename}"
        
        if [[ "${SKIP_EXISTING}" -eq 1 && -f "${filepath}" ]]; then
            local actual_size
            actual_size=$(get_file_size "${filepath}")
            
            if [[ "${actual_size}" -eq "${expected_size}" ]]; then
                ((skipped++))
                continue
            fi
        fi
        
        echo "${url}" >> "00_meta/${ACCESSION}.urls_to_download.txt"
        
    done < "00_meta/${ACCESSION}.urls.md5.bytes.tsv"
    
    local to_download
    to_download=$(wc -l < "00_meta/${ACCESSION}.urls_to_download.txt" | tr -d ' ')
    
    if [[ ${to_download} -eq 0 ]]; then
        log_info "All ${file_count} files already exist with correct size"
    else
        log_info "Downloading ${to_download} files (${skipped} existing skipped)"
        log_info "Using ${PARALLEL_DL} parallel connections..."
        
        set +e
        xargs -n 1 -P "${PARALLEL_DL}" -I {} \
            wget -c -nv \
                 -t "${WGET_RETRIES}" \
                 --waitretry="${WGET_WAIT}" \
                 --timeout="${WGET_TIMEOUT}" \
                 --read-timeout="${WGET_READ_TIMEOUT}" \
                 -P "01_fastq" \
                 "{}" \
            < "00_meta/${ACCESSION}.urls_to_download.txt" 2>&1 | tee "99_logs/wget.log"
        set -e
    fi
    
    local dl_end
    dl_end=$(date +%s)
    local dl_duration=$((dl_end - dl_start))
    
    log_info "Download completed in $(human_time ${dl_duration})"
    
    #---------------------------------------------------------------------------
    # Step 5: 検証
    #---------------------------------------------------------------------------
    log_header "Step 5: Verifying Files"
    
    local verify_start
    verify_start=$(date +%s)
    
    local current=0
    local verified=0
    local failed=0
    local failed_files=()
    
    while IFS=$'\t' read -r url expected_md5 expected_bytes run_acc; do
        ((current++))
        
        local filename
        filename="$(basename "${url}")"
        local filepath="01_fastq/${filename}"
        
        printf "\r[%d/%d] %s" "${current}" "${file_count}" "${filename:0:40}"
        
        # 存在確認
        if [[ ! -f "${filepath}" ]]; then
            failed_files+=("${filename}")
            ((failed++))
            continue
        fi
        
        # サイズ確認
        local actual_size
        actual_size=$(get_file_size "${filepath}")
        
        if [[ "${actual_size}" -eq 0 ]]; then
            failed_files+=("${filename}")
            ((failed++))
            continue
        fi
        
        if [[ -n "${expected_bytes}" && "${actual_size}" -ne "${expected_bytes}" ]]; then
            failed_files+=("${filename}")
            ((failed++))
            continue
        fi
        
        # MD5検証
        if [[ "${VERIFY_MD5}" -eq 1 && -n "${expected_md5}" ]]; then
            local actual_md5
            actual_md5=$(calculate_md5 "${filepath}")
            
            if [[ "${actual_md5}" != "${expected_md5}" ]]; then
                failed_files+=("${filename}")
                ((failed++))
                continue
            fi
        fi
        
        ((verified++))
        
    done < "00_meta/${ACCESSION}.urls.md5.bytes.tsv"
    
    printf "\r%-60s\n" ""
    
    local verify_end
    verify_end=$(date +%s)
    local verify_duration=$((verify_end - verify_start))
    
    local verify_type="SIZE"
    [[ "${VERIFY_MD5}" -eq 1 ]] && verify_type="MD5"
    
    log_info "Verification (${verify_type}) completed in $(human_time ${verify_duration})"
    log_info "Verified: ${verified}/${file_count}"
    
    if [[ ${failed} -gt 0 ]]; then
        log_error "Failed: ${failed}/${file_count}"
    fi
    
    #---------------------------------------------------------------------------
    # 結果サマリー
    #---------------------------------------------------------------------------
    local script_end
    script_end=$(date +%s)
    local total_duration=$((script_end - script_start))
    
    log_header "Summary"
    
    if [[ ${failed} -gt 0 ]]; then
        log_error "INCOMPLETE: ${failed} files failed verification"
        echo ""
        for f in "${failed_files[@]}"; do
            log_error "  - ${f}"
        done
        
        printf '%s\n' "${failed_files[@]}" > "99_logs/failed_files.txt"
        
        echo ""
        log_info "Failed files saved to: 99_logs/failed_files.txt"
        log_info ""
        log_info "To retry:"
        log_info "  cd ${work_dir}"
        log_info "  for f in \$(cat 99_logs/failed_files.txt); do"
        log_info "    grep \"\${f}\" 00_meta/${ACCESSION}.urls.txt | xargs wget -c -P 01_fastq"
        log_info "  done"
        
        exit 1
    fi
    
    log_success "All ${verified} files verified successfully!"
    echo ""
    log_info "Accession:   ${ACCESSION}"
    log_info "Location:    ${work_dir}/01_fastq"
    log_info "Files:       ${verified}"
    log_info "Total size:  $(human_bytes "${total_bytes}")"
    log_info "Verified:    ${verify_type}"
    echo ""
    log_info "Time: Download $(human_time ${dl_duration}), Verify $(human_time ${verify_duration}), Total $(human_time ${total_duration})"
    
    # 完了マーカー
    {
        echo "completed: $(date '+%Y-%m-%d %H:%M:%S')"
        echo "accession: ${ACCESSION}"
        echo "files: ${verified}"
        echo "bytes: ${total_bytes}"
        echo "md5_verified: ${VERIFY_MD5}"
        echo "duration_seconds: ${total_duration}"
    } > "99_logs/completed.txt"
    
    echo ""
    log_success "DONE"
}

#-------------------------------------------------------------------------------
# エントリーポイント
#-------------------------------------------------------------------------------
main "$@"