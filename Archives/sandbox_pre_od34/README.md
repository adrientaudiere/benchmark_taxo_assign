Sandbox scripts that read targets from before OD 34 (`d_all_taxo`,
`d_asv_for_assignation`, `d_asv_common`): the raw-ASV inputs, which the
stores no longer build and which `tar_prune()` removes. Moved here on
2026-10-07 (code review §1.8) rather than ported, so that the figures
`HISTORY.md` cites from them keep the meaning they had when measured.
They no longer run as they are (their `here::i_am()` names the old path).
