# =============================================================================
# _genotyping_configs.sh — shared ambientmapper genotyping config helper
#
# This is a SOURCED bash helper. It is NOT meant to be executed directly.
# Wrappers source it to obtain the C0 baseline flag block, the case statement
# of all 55 supported genotyping configurations (Phase 1's 42 C0..C4d configs
# plus Phase 3's 13 S01..S13 stacked configs), the friend-rescue routing
# logic, and the ambientmapper invocation wrapper.
#
# WHO USES THIS:
#   - workflows/03_genotyping/synthetic/03_21_genotyping_synthetic_factorial.sh   (Phase 2, Fig. 4A to C, S5)
#   - workflows/03_genotyping/root1_sub1k/02_11_genotyping_stacked_phase3_sub1k.sh  (Phase 3, Fig. S6)
#   - workflows/03_genotyping/root1/02_12_genotyping_full_root1_phase4.sh     (Phase 4, Fig. 4D to G)
#   - workflows/03_genotyping/sm2v2/01_09_genotyping_SM2v2.sh                  (SM2v2 C0 run, Fig. 2, 3, 5)
#   - workflows/03_genotyping/zhang2024/04_05a_genotyping_4cfg_B73Mo17.sh and
#     04_05b_genotyping_4cfg_multi.sh                                          (Zhang 2024 C0 runs, Fig. 4H to N)
#
# WHO DOES NOT:
#   - workflows/03_genotyping/root1_sub1k/02_10_genotyping_factorial_sub1k.sh — Phase 1, kept inline
#     to preserve the ability to re-run any Phase 1 task without diff risk.
#     Accepting the duplication is intentional.
#
# CONVENTION:
#   The leading underscore marks it as "library, not entry point".
#
# WHY A FUNCTION-BASED HELPER (and not a vars-only file):
#   apply_config_overrides expects to mutate caller-scope variables. Bash has
#   no module system; sourcing into the caller scope is the standard pattern.
#
# CRITICAL — bash word-splitting:
#   The wrapper run_ambientmapper_genotyping uses UNQUOTED expansion for
#   ${XA_MAX_FLAG}, ${RECLASS_FLAG}, ${WDISC_FLAG}, ${WDISC_MODE_FLAG},
#   ${XMAP_FLAG}. These variables can be the empty string (e.g. C3c_xa_unlim
#   sets XA_MAX_FLAG=""). Unquoted expansion makes the empty case expand to
#   ZERO tokens; quoting them would pass a literal empty string to
#   ambientmapper and break the call. Test the empty case before promoting
#   any change to this file.
# =============================================================================

# Guard against direct execution.
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "ERROR: _genotyping_configs.sh must be sourced, not executed." >&2
    echo "       Use:  source workflows/03_genotyping/_genotyping_configs.sh" >&2
    exit 1
fi

# -----------------------------------------------------------------------------
# init_c0_defaults
#
# Sets the C0 baseline scalars in the caller's scope. Mirrors lines 191-224
# of 02_10_genotyping_factorial_sub1k.sh exactly. Call this once at the top
# of every wrapper, BEFORE apply_config_overrides.
# -----------------------------------------------------------------------------
init_c0_defaults() {
    # Read-quality filters
    MAPQ_MIN=20
    XA_MAX_FLAG="--max-xa-genotyping 0"     # empty string disables the XA filter

    # Topk + scoring weights
    TOPK=2
    W_CONFIDENT=1.0
    W_AMBIG=0.1
    W_AS=1.0
    W_MAPQ=1.0
    W_NM=1.0

    # BIC thresholds (singlet/doublet)
    MASS_MIN=0.6
    RATIO_MIN=1.5
    MAX_ALPHA=0.5
    BIC_MARGIN=6
    DOUBLET_MINOR_MIN=0.20

    # Eta (ambient profile)
    ETA_ITERS=2

    # Empty-cell thresholds (frozen across the factorial)
    EMPTY_BIC_MARGIN=6
    EMPTY_TOP1_MAX=0.6
    EMPTY_RATIO12_MAX=2
    EMPTY_SEED_BIC_MIN=10
    EMPTY_TAU_QUANTILE=0.95

    # Rescue features (default = ALL ON, except xmap)
    RECLASS_FLAG="--topk-reclass"
    WDISC_FLAG="--winner-discount"
    WDISC_MODE_FLAG="--winner-discount-mode winner_ratio"
    XMAP_FLAG="--no-xmap"
}

# -----------------------------------------------------------------------------
# apply_config_overrides <CONFIG_NAME>
#
# Overrides the C0 defaults set by init_c0_defaults. Each config branch only
# sets what differs from C0; everything else stays at the baseline value.
#
# Supported configs (55):
#   - 42 Phase 1: C0, C1a..C1g, C2a..C2e, C3a..C3e, C4a* (10), C4b* (6),
#     C4c_eta0/eta5, C4d_wamb*_wd{on,off} (6)
#   - 13 Phase 3: S01..S13 stacked configurations
#
# Exits on unknown config.
# -----------------------------------------------------------------------------
apply_config_overrides() {
    local CONFIG="$1"
    if [[ -z "${CONFIG}" ]]; then
        echo "ERROR: apply_config_overrides requires a CONFIG name argument" >&2
        return 1
    fi

    case "${CONFIG}" in
        # ---------------------------------------------------------------------
        # Phase 1 C0 — pure baseline
        # ---------------------------------------------------------------------
        C0)
            : # no overrides
            ;;

        # ---------------------------------------------------------------------
        # Phase 1 C1 — rescue knockout factorial
        # ---------------------------------------------------------------------
        C1a_noreclass)
            RECLASS_FLAG="--no-topk-reclass"
            ;;
        C1b_nowdisc)
            WDISC_FLAG="--no-winner-discount"
            ;;
        C1c_nofriend)
            : # friend rescue knocked out via the friendwithout chunks dir
            ;;
        C1d_noreclass_nowdisc)
            RECLASS_FLAG="--no-topk-reclass"
            WDISC_FLAG="--no-winner-discount"
            ;;
        C1e_noreclass_nofriend)
            RECLASS_FLAG="--no-topk-reclass"
            ;;
        C1f_nowdisc_nofriend)
            WDISC_FLAG="--no-winner-discount"
            ;;
        C1g_naked)
            RECLASS_FLAG="--no-topk-reclass"
            WDISC_FLAG="--no-winner-discount"
            ;;

        # ---------------------------------------------------------------------
        # Phase 1 C2 — xmap layered + fixeta-ma08 disentanglement
        # ---------------------------------------------------------------------
        C2a_xmap)
            XMAP_FLAG=""   # omit --no-xmap → xmap enabled
            ;;
        C2b_xmap_noreclass)
            XMAP_FLAG=""
            RECLASS_FLAG="--no-topk-reclass"
            ;;
        C2c_xmap_eta0)
            XMAP_FLAG=""
            ETA_ITERS=0
            ;;
        C2d_xmap_ma08)
            XMAP_FLAG=""
            MAX_ALPHA=0.8
            ;;
        C2e_xmap_eta0_ma08)
            XMAP_FLAG=""
            ETA_ITERS=0
            MAX_ALPHA=0.8
            ;;

        # ---------------------------------------------------------------------
        # Phase 1 C3 — read filter ablations
        # ---------------------------------------------------------------------
        C3a_mq10)
            MAPQ_MIN=10
            ;;
        C3b_mq50)
            MAPQ_MIN=50
            ;;
        C3c_xa_unlim)
            XA_MAX_FLAG=""   # omit the XA filter entirely
            ;;
        C3d_mq10_xa_unlim)
            MAPQ_MIN=10
            XA_MAX_FLAG=""
            ;;
        C3e_mq50_xa_unlim)
            MAPQ_MIN=50
            XA_MAX_FLAG=""
            ;;

        # ---------------------------------------------------------------------
        # Ad-hoc: C0 + xmap + mq50 (off-factorial 2x2 cell completing the
        # {xmap on/off} × {mq50 on/off} grid alongside C0, C2a_xmap, C3b_mq50).
        # Used by 04_05a_genotyping_4cfg_B73Mo17.sh and 04_05b_genotyping_4cfg_multi.sh
        # (only the C0 cell of that grid is displayed in the manuscript).
        # Friend rescue stays ON (default), reclass/wdisc stay ON, eta_iters
        # stays at C0 default = 2 (matches C2a_xmap, NOT S03/S04 which set eta=0).
        # ---------------------------------------------------------------------
        Cxmap_mq50)
            MAPQ_MIN=50
            XMAP_FLAG=""   # omit --no-xmap → xmap enabled
            ;;

        # ---------------------------------------------------------------------
        # Phase 1 C4a — BIC singlet thresholds
        # ---------------------------------------------------------------------
        C4a_mass04)  MASS_MIN=0.4 ;;
        C4a_mass05)  MASS_MIN=0.5 ;;
        C4a_mass07)  MASS_MIN=0.7 ;;
        C4a_ratio12) RATIO_MIN=1.2 ;;
        C4a_ratio18) RATIO_MIN=1.8 ;;
        C4a_ratio20) RATIO_MIN=2.0 ;;
        C4a_ma03)    MAX_ALPHA=0.3 ;;
        C4a_ma08)    MAX_ALPHA=0.8 ;;
        C4a_bic3)    BIC_MARGIN=3 ;;
        C4a_bic10)   BIC_MARGIN=10 ;;

        # ---------------------------------------------------------------------
        # Phase 1 C4b — scoring weights
        # ---------------------------------------------------------------------
        C4b_was05) W_AS=0.5 ;;
        C4b_was15) W_AS=1.5 ;;
        C4b_wmq05) W_MAPQ=0.5 ;;
        C4b_wmq20) W_MAPQ=2.0 ;;
        C4b_wnm05) W_NM=0.5 ;;
        C4b_wnm20) W_NM=2.0 ;;

        # ---------------------------------------------------------------------
        # Phase 1 C4c — eta iters
        # ---------------------------------------------------------------------
        C4c_eta0) ETA_ITERS=0 ;;
        C4c_eta5) ETA_ITERS=5 ;;

        # ---------------------------------------------------------------------
        # Phase 1 C4d — w_ambiguous × wdisc structural cross
        # ---------------------------------------------------------------------
        C4d_wamb005_wdon)  W_AMBIG=0.05 ;;
        C4d_wamb005_wdoff) W_AMBIG=0.05; WDISC_FLAG="--no-winner-discount" ;;
        C4d_wamb02_wdon)   W_AMBIG=0.2 ;;
        C4d_wamb02_wdoff)  W_AMBIG=0.2;  WDISC_FLAG="--no-winner-discount" ;;
        C4d_wamb05_wdon)   W_AMBIG=0.5 ;;
        C4d_wamb05_wdoff)  W_AMBIG=0.5;  WDISC_FLAG="--no-winner-discount" ;;

        # ---------------------------------------------------------------------
        # Phase 3 — stacked configurations (S01..S13)
        #
        # Reference table (mq, xa, friend, xmap, eta, w_amb, bic_margin):
        #   S01: 50,  0, off, OFF, 2, 0.1, 6   minimal stack (mq50 + nofriend)
        #   S02: 50,  ∞, off, OFF, 2, 0.1, 6   + XA unlimited
        #   S03: 50,  0, off, ON,  0, 0.1, 6   + xmap (with eta0)
        #   S04: 50,  ∞, off, ON,  0, 0.1, 6   + xmap + XA unlim
        #   S05: 50,  0, off, OFF, 2, 0.5, 6   + w_amb=0.5
        #   S06: 50,  ∞, off, OFF, 2, 0.5, 6   + w_amb=0.5 + XA unlim
        #   S07: 50,  0, off, ON,  0, 0.5, 6   + xmap + w_amb=0.5
        #   S08: 50,  0, off, ON,  0, 0.5, 6   FULL STACK with XA=0
        #   S09: 50,  ∞, off, ON,  0, 0.5, 6   FULL STACK with XA unlim
        #   S10: 50,  0, off, ON,  0, 0.5, 3   FULL STACK + bic_margin=3 bonus
        #   S11: 50,  0, ON,  ON,  0, 0.5, 6   sanity: friend ON inside stack
        #   S12: 50,  0, off, ON,  0, 0.1, 6   sanity: w_amb=0.1 inside stack
        #   S13: 50,  0, off, ON,  2, 0.5, 6   sanity: eta_iters=2 inside stack
        #
        # Note: friend OFF is implemented via the friendwithout chunks dir,
        # NOT via a flag. friend_mode_for_config encodes the routing.
        # ---------------------------------------------------------------------
        S01_nofr_mq50)
            MAPQ_MIN=50
            ;;
        S02_nofr_mq50_xaUL)
            MAPQ_MIN=50
            XA_MAX_FLAG=""
            ;;
        S03_nofr_mq50_xmap)
            MAPQ_MIN=50
            XMAP_FLAG=""
            ETA_ITERS=0
            ;;
        S04_nofr_mq50_xmap_xaUL)
            MAPQ_MIN=50
            XA_MAX_FLAG=""
            XMAP_FLAG=""
            ETA_ITERS=0
            ;;
        S05_nofr_mq50_wamb)
            MAPQ_MIN=50
            W_AMBIG=0.5
            ;;
        S06_nofr_mq50_wamb_xaUL)
            MAPQ_MIN=50
            W_AMBIG=0.5
            XA_MAX_FLAG=""
            ;;
        S07_nofr_mq50_xmap_wamb)
            MAPQ_MIN=50
            XMAP_FLAG=""
            ETA_ITERS=0
            W_AMBIG=0.5
            ;;
        S08_full_xa0)
            MAPQ_MIN=50
            XMAP_FLAG=""
            ETA_ITERS=0
            W_AMBIG=0.5
            ;;
        S09_full_xaUL)
            MAPQ_MIN=50
            XA_MAX_FLAG=""
            XMAP_FLAG=""
            ETA_ITERS=0
            W_AMBIG=0.5
            ;;
        S10_full_bic3)
            MAPQ_MIN=50
            XMAP_FLAG=""
            ETA_ITERS=0
            W_AMBIG=0.5
            BIC_MARGIN=3
            ;;
        S11_full_friendon)
            # NOTE: friend stays ON via friend_mode_for_config; this case
            # only sets the genotyping flags. Friend mode is determined by
            # the chunks dir routing (friend=with).
            MAPQ_MIN=50
            XMAP_FLAG=""
            ETA_ITERS=0
            W_AMBIG=0.5
            ;;
        S12_full_wamb01)
            # Sanity: full stack but w_amb stays at 0.1 (the C0 default)
            MAPQ_MIN=50
            XMAP_FLAG=""
            ETA_ITERS=0
            # W_AMBIG stays at C0 default 0.1
            ;;
        S13_full_eta2)
            # Sanity: full stack but eta_iters stays at 2 (the C0 default)
            MAPQ_MIN=50
            XMAP_FLAG=""
            # ETA_ITERS stays at C0 default 2
            W_AMBIG=0.5
            ;;

        *)
            echo "ERROR: unknown config '${CONFIG}'" >&2
            echo "       see _genotyping_configs.sh for the list of supported names" >&2
            return 1
            ;;
    esac
}

# -----------------------------------------------------------------------------
# friend_mode_for_config <CONFIG_NAME>
#
# Echoes 'with' or 'without' depending on whether the config requires the
# friend-rescue-knocked-out chunks dir.
#
# Phase 1: C1c, C1e, C1f, C1g use 'without'.
# Phase 3: S01-S10, S12, S13 use 'without'. Only S11 uses 'with'.
#
# Usage:
#   FRIEND=$(friend_mode_for_config "$CONFIG")
# -----------------------------------------------------------------------------
friend_mode_for_config() {
    local CONFIG="$1"
    case "${CONFIG}" in
        # Phase 1 friend-knockout configs
        C1c_nofriend|C1e_noreclass_nofriend|C1f_nowdisc_nofriend|C1g_naked)
            echo "without"
            ;;
        # Phase 3 stacked configs — explicitly listed (no fall-through)
        S01_nofr_mq50|S02_nofr_mq50_xaUL|S03_nofr_mq50_xmap|S04_nofr_mq50_xmap_xaUL)
            echo "without"
            ;;
        S05_nofr_mq50_wamb|S06_nofr_mq50_wamb_xaUL|S07_nofr_mq50_xmap_wamb)
            echo "without"
            ;;
        S08_full_xa0|S09_full_xaUL|S10_full_bic3)
            echo "without"
            ;;
        S11_full_friendon)
            echo "with"
            ;;
        S12_full_wamb01|S13_full_eta2)
            echo "without"
            ;;
        # All Phase 1 non-knockout configs
        *)
            echo "with"
            ;;
    esac
}

# -----------------------------------------------------------------------------
# print_effective_flags
#
# Prints the current C0+overrides flag block to stdout. Call this AFTER
# init_c0_defaults && apply_config_overrides "$CONFIG" to log the effective
# flags for a run. Mirrors lines 369-381 of 02_10.
# -----------------------------------------------------------------------------
print_effective_flags() {
    echo "  effective flags:"
    echo "    --min-mapq-genotyping ${MAPQ_MIN}"
    echo "    ${XA_MAX_FLAG:-(no XA filter)}"
    echo "    --topk-genomes ${TOPK}"
    echo "    --w-confident ${W_CONFIDENT} --w-ambiguous ${W_AMBIG}"
    echo "    --w-as ${W_AS} --w-mapq ${W_MAPQ} --w-nm ${W_NM}"
    echo "    --single-mass-min ${MASS_MIN} --ratio-top1-top2-min ${RATIO_MIN}"
    echo "    --max-alpha ${MAX_ALPHA} --bic-margin ${BIC_MARGIN}"
    echo "    --eta-iters ${ETA_ITERS}"
    echo "    ${RECLASS_FLAG}  ${WDISC_FLAG}  ${WDISC_MODE_FLAG}  ${XMAP_FLAG}"
}

# -----------------------------------------------------------------------------
# check_ambientmapper_cli_features
#
# Verifies the installed ambientmapper has --winner-discount and --topk-reclass
# flags via Typer/Click introspection. Mirrors 02_10:347-360. Use this in
# wrappers before calling run_ambientmapper_genotyping.
# -----------------------------------------------------------------------------
check_ambientmapper_cli_features() {
    local _GENO_FLAGS
    _GENO_FLAGS=$(python -c "import ambientmapper.cli as _c, typer; g = typer.main.get_command(_c.app).commands['genotyping']; names = set(); [names.update(p.opts) for p in g.params]; print(' '.join(sorted(names)))" 2>/dev/null)
    if [[ -z "${_GENO_FLAGS}" ]]; then
        echo "ERROR: failed to introspect ambientmapper.cli.app — is the package installed?" >&2
        return 1
    fi
    if [[ ! " ${_GENO_FLAGS} " == *" --winner-discount "* ]]; then
        echo "ERROR: --winner-discount flag not found in installed ambientmapper" >&2
        return 1
    fi
    if [[ ! " ${_GENO_FLAGS} " == *" --topk-reclass "* ]]; then
        echo "ERROR: --topk-reclass flag not found in installed ambientmapper" >&2
        return 1
    fi
}

# -----------------------------------------------------------------------------
# run_ambientmapper_genotyping <CONFIG_JSON> <ASSIGN_GLOB> <OUTDIR>
#
# Runs ambientmapper genotyping with the currently-set flag block. Mirrors
# the invocation at 02_10:384-417 EXACTLY. The empty-string flag variables
# (XA_MAX_FLAG, RECLASS_FLAG, WDISC_FLAG, WDISC_MODE_FLAG, XMAP_FLAG) are
# expanded UNQUOTED to allow them to expand to zero tokens when empty.
#
# Caller must have already called:
#   init_c0_defaults
#   apply_config_overrides "$CONFIG"
#
# and must have created OUTDIR.
# -----------------------------------------------------------------------------
run_ambientmapper_genotyping() {
    local CONFIG_JSON="$1"
    local ASSIGN_GLOB="$2"
    local OUTDIR="$3"

    if [[ -z "${CONFIG_JSON}" || -z "${ASSIGN_GLOB}" || -z "${OUTDIR}" ]]; then
        echo "ERROR: run_ambientmapper_genotyping needs 3 args: CONFIG_JSON ASSIGN_GLOB OUTDIR" >&2
        return 1
    fi
    if [[ ! -f "${CONFIG_JSON}" ]]; then
        echo "ERROR: config JSON not found: ${CONFIG_JSON}" >&2
        return 1
    fi

    # Winner-only mode is OFF (--no-winner-only) and beta = 10 in every run that
    # goes through this helper. The library defaults are winner-only ON and
    # beta = 1. Neither flag is echoed by print_effective_flags, so the two lines
    # below are the only record of the values used.
    # shellcheck disable=SC2086
    ambientmapper genotyping \
        --config "${CONFIG_JSON}" \
        --assign "${ASSIGN_GLOB}" \
        --outdir "${OUTDIR}" \
        --no-resume \
        --threads 8 \
        --pass1-workers 8 \
        --no-winner-only \
        --beta 10 \
        --min-reads 5 \
        --chunk-rows 100000 \
        --w-confident ${W_CONFIDENT} \
        --w-ambiguous ${W_AMBIG} \
        --w-as ${W_AS} \
        --w-mapq ${W_MAPQ} \
        --w-nm ${W_NM} \
        --topk-genomes ${TOPK} \
        --doublet-minor-min ${DOUBLET_MINOR_MIN} \
        --single-mass-min ${MASS_MIN} \
        --ratio-top1-top2-min ${RATIO_MIN} \
        --max-alpha ${MAX_ALPHA} \
        --bic-margin ${BIC_MARGIN} \
        --eta-iters ${ETA_ITERS} \
        --empty-bic-margin ${EMPTY_BIC_MARGIN} \
        --empty-top1-max ${EMPTY_TOP1_MAX} \
        --empty-ratio12-max ${EMPTY_RATIO12_MAX} \
        --empty-seed-bic-min ${EMPTY_SEED_BIC_MIN} \
        --empty-tau-quantile ${EMPTY_TAU_QUANTILE} \
        --min-mapq-genotyping ${MAPQ_MIN} \
        ${XA_MAX_FLAG} \
        ${RECLASS_FLAG} \
        ${WDISC_FLAG} \
        ${WDISC_MODE_FLAG} \
        ${XMAP_FLAG}
}
