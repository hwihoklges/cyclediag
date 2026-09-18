"""Detect RPT / capacheck cycles and build routine-life exclusion sets."""



from __future__ import annotations



from dataclasses import dataclass, field



import numpy as np

import pandas as pd



# Routine SJ900 life cycle: charge – rest – discharge – rest (4 step-end rows)

ROUTINE_STEPS = 4

POST_RPT_EXCLUDE = 5

DCHG_ROUTINE_MIN_AH = 50.0

# Full 0.33C capa vs 0.5C routine: |I_capa| ≈ 0.67 · |I_routine| (e.g. 25.8 vs 38.7 A)

CAPA_FULL_I_RATIO = 0.85

CAPA_FULL_MIN_AH = 10.0





@dataclass

class ProtocolExclusion:

    """Cycles excluded from routine peak / fade analysis."""



    excluded: set[int] = field(default_factory=set)

    rpt_cycles: set[int] = field(default_factory=set)

    capacheck_cycles: set[int] = field(default_factory=set)

    capa_full_cycles: set[int] = field(default_factory=set)

    post_rpt_cycles: set[int] = field(default_factory=set)

    rpt_blocks: list[list[int]] = field(default_factory=list)

    flags: pd.DataFrame = field(default_factory=pd.DataFrame)

    post_rpt_exclude: int = POST_RPT_EXCLUDE



    def to_meta(self) -> dict:

        return {

            "routine_steps": ROUTINE_STEPS,

            "post_rpt_exclude": self.post_rpt_exclude,

            "n_rpt": len(self.rpt_cycles),

            "n_capacheck": len(self.capacheck_cycles),

            "n_capa_full": len(self.capa_full_cycles),

            "n_post_rpt": len(self.post_rpt_cycles),

            "n_excluded_total": len(self.excluded),

            "rpt_blocks": self.rpt_blocks,

            "excluded_cycles": sorted(self.excluded),

        }





def _dchg_ah(grp: pd.DataFrame) -> float:

    caps = pd.to_numeric(grp.get("discharge_capacity"), errors="coerce")

    if caps is None or not caps.notna().any():

        return 0.0

    # Normalization owns raw-unit conversion; bare logical columns are Ah.
    return float(caps.max())





def _discharge_abs_i(grp: pd.DataFrame) -> float:

    """Median |I| on discharge steps (ignores rare outliers in StepEnd)."""

    cur = pd.to_numeric(grp.get("current"), errors="coerce").fillna(0.0)

    disc = cur[cur < -1.0]

    if disc.empty:

        return 0.0

    return float(disc.abs().median())





def detect_protocol_flags(

    step_df: pd.DataFrame,

    *,

    routine_steps: int = ROUTINE_STEPS,

    dchg_routine_min_ah: float = DCHG_ROUTINE_MIN_AH,

    capa_full_i_ratio: float = CAPA_FULL_I_RATIO,

    capa_full_min_ah: float = CAPA_FULL_MIN_AH,

) -> pd.DataFrame:

    """Per-cycle flags from step-end table (one row per protocol step)."""

    if step_df.empty or "cycle" not in step_df.columns:

        return pd.DataFrame()



    rows: list[dict] = []

    for cyc, grp in step_df.groupby("cycle"):

        cyc_i = int(cyc)

        n_steps = len(grp)

        dchg_ah = _dchg_ah(grp)

        cur = pd.to_numeric(grp.get("current"), errors="coerce").fillna(0.0)

        i_unique = int(cur.round(3).nunique())

        i_dchg_abs = _discharge_abs_i(grp)



        is_non_routine_steps = n_steps != routine_steps

        is_low_dchg = 0 < dchg_ah < dchg_routine_min_ah

        # Partial SOC / DC-IR steps inside RPT block

        is_capacheck = is_non_routine_steps and is_low_dchg

        is_rpt = is_non_routine_steps



        rows.append({

            "cycle": cyc_i,

            "n_steps": n_steps,

            "dchg_Ah": round(dchg_ah, 4),

            "i_dchg_abs": round(i_dchg_abs, 4),

            "i_unique": i_unique,

            "is_non_routine_steps": is_non_routine_steps,

            "is_capacheck": is_capacheck,

            "is_capa_full": False,  # filled below from rate vs routine

            "is_rpt": is_rpt,

            "protocol_kind": (

                "rpt+capacheck" if is_non_routine_steps and is_capacheck

                else "rpt" if is_non_routine_steps

                else "capacheck" if is_capacheck

                else "routine"

            ),

        })



    flags = pd.DataFrame(rows).sort_values("cycle").reset_index(drop=True)

    if flags.empty:

        return flags



    # Modal discharge |I| among 4-step, full-capacity candidates = 0.5C routine rate

    four = flags[

        (flags["n_steps"] == routine_steps)

        & (flags["dchg_Ah"] >= dchg_routine_min_ah)

        & (flags["i_dchg_abs"] > 1.0)

    ]

    if len(four) >= 3:

        # Round to 0.5 A bins so 38.67 / 38.67 cluster; pick the most common bin center

        bins = (four["i_dchg_abs"] / 0.5).round() * 0.5

        mode_i = float(bins.mode().iloc[0])

        capa_mask = (

            (flags["n_steps"] == routine_steps)

            & (flags["dchg_Ah"] >= capa_full_min_ah)

            & (flags["i_dchg_abs"] > 1.0)

            & (flags["i_dchg_abs"] < mode_i * capa_full_i_ratio)

        )

        flags.loc[capa_mask, "is_capa_full"] = True

        flags.loc[capa_mask, "is_capacheck"] = True

        flags.loc[capa_mask, "protocol_kind"] = "capa_full"



    return flags





def rpt_blocks(sorted_rpt: list[int]) -> list[list[int]]:

    if not sorted_rpt:

        return []

    blocks: list[list[int]] = []

    blk = [sorted_rpt[0]]

    for c in sorted_rpt[1:]:

        if c == blk[-1] + 1:

            blk.append(c)

        else:

            blocks.append(blk)

            blk = [c]

    blocks.append(blk)

    return blocks





def preceding_capa_full_cycles(

    flags: pd.DataFrame,

    block_start: int,

) -> list[int]:

    """Consecutive 0.33C full-capa cycles immediately before an RPT/DC-IR block."""

    if flags.empty or "is_capa_full" not in flags.columns:

        return []

    by_cyc = flags.set_index("cycle")

    out: list[int] = []

    c = int(block_start) - 1

    while c in by_cyc.index and bool(by_cyc.loc[c, "is_capa_full"]):

        out.append(c)

        c -= 1

    return sorted(out)





def build_protocol_exclusion(

    step_df: pd.DataFrame,

    *,

    post_rpt_exclude: int = POST_RPT_EXCLUDE,

    routine_steps: int = ROUTINE_STEPS,

    dchg_routine_min_ah: float = DCHG_ROUTINE_MIN_AH,

) -> ProtocolExclusion:

    """RPT/capacheck cycles + N recovery cycles after each RPT block."""

    flags = detect_protocol_flags(

        step_df,

        routine_steps=routine_steps,

        dchg_routine_min_ah=dchg_routine_min_ah,

    )

    if flags.empty:

        return ProtocolExclusion(flags=flags, post_rpt_exclude=post_rpt_exclude)



    rpt_cycles = set(flags.loc[flags["is_rpt"], "cycle"].astype(int))

    capa_full_cycles = set(flags.loc[flags["is_capa_full"], "cycle"].astype(int))

    capacheck_cycles = set(flags.loc[flags["is_capacheck"], "cycle"].astype(int))

    # Include full 0.33C capa in capacheck set for checkpoint preference / reporting

    capacheck_cycles |= capa_full_cycles

    blocks = rpt_blocks(sorted(rpt_cycles))



    excluded = set(rpt_cycles) | capa_full_cycles

    post_rpt: set[int] = set()

    for block in blocks:

        end = block[-1]

        for k in range(1, post_rpt_exclude + 1):

            post_rpt.add(end + k)

    excluded |= post_rpt



    flags = flags.copy()

    flags["protocol_excluded"] = flags["cycle"].isin(excluded)

    flags["post_rpt_buffer"] = flags["cycle"].isin(post_rpt)

    only_post = flags["post_rpt_buffer"] & ~flags["is_rpt"] & ~flags["is_capa_full"]

    flags.loc[only_post, "protocol_kind"] = "post_rpt"



    return ProtocolExclusion(

        excluded=excluded,

        rpt_cycles=rpt_cycles,

        capacheck_cycles=capacheck_cycles,

        capa_full_cycles=capa_full_cycles,

        post_rpt_cycles=post_rpt,

        rpt_blocks=blocks,

        flags=flags,

        post_rpt_exclude=post_rpt_exclude,

    )





def apply_protocol_exclusion(

    wide_df: pd.DataFrame,

    long_df: pd.DataFrame,

    protocol: ProtocolExclusion,

) -> tuple[pd.DataFrame, pd.DataFrame]:

    """Mark RPT/capacheck cycles non-usable on wide/long peak tables."""

    if wide_df.empty or not protocol.excluded:

        return wide_df, long_df



    kind_map: dict[int, str] = {}

    if not protocol.flags.empty:

        kind_map = dict(zip(protocol.flags["cycle"].astype(int), protocol.flags["protocol_kind"]))



    wide = wide_df.copy()

    wide["protocol_excluded"] = wide["cycle"].isin(protocol.excluded)

    wide["protocol_kind"] = wide["cycle"].map(kind_map).fillna("routine")



    base_flags = wide["exclude_flags"].fillna("").astype(str)

    tags = "protocol:" + wide["protocol_kind"].astype(str)

    wide["exclude_flags"] = np.where(

        wide["protocol_excluded"],

        np.where(base_flags == "", tags, base_flags + "|" + tags),

        base_flags,

    )

    wide["usable_auto"] = wide["usable_auto"] & ~wide["protocol_excluded"]

    wide["usable"] = wide["usable_auto"]



    usable_map = dict(zip(wide["cycle"], wide["usable"]))

    cha_map = dict(zip(wide["cycle"], wide["usable_charge"]))

    dis_map = dict(zip(wide["cycle"], wide["usable_discharge"]))



    long = long_df.copy()

    long["protocol_excluded"] = long["cycle"].isin(protocol.excluded)

    long["protocol_kind"] = long["cycle"].map(kind_map).fillna("routine")

    long["usable"] = long["cycle"].map(usable_map)

    is_charge = long["leg"] == "charge"

    long["usable_leg"] = np.where(is_charge, long["cycle"].map(cha_map), long["cycle"].map(dis_map))

    return wide, long


