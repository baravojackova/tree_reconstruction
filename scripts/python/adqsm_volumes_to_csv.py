# -*- coding: utf-8 -*-
# =====================================================================
#  Zapíše objemy AdQSM (z TreesParams.txt) do tabulky výsledků dané série.
# ---------------------------------------------------------------------
#  K čemu to je: adtree_reconstruct_compare.py umí zapsat řádky AdQSM jen
#  spolu s AdTree (potřebuje skeleton .ply). U série B24 AdTree zatím
#  neděláme, proto tenhle malý skript dělá JEN AdQSM část.
#
#  Co skript dělá, pro každý strom a každou variantu:
#   1) najde data/<série>/<strom>/<varianta>/TreesParams.txt,
#   2) přečte z něj objemy, DBH, výšku a délky (parse_adqsm_params_file),
#   3) zapíše jeden řádek do volume_results_<série>.csv (upsert_result).
#  Řádek se pozná podle dvojice (strom, název metody). Opakované spuštění
#  řádek PŘEPÍŠE, nevznikne duplicita.
#
#  Řádky mají PŘESNĚ stejný název metody a stejné sloupce jako ty, které
#  zapisuje adtree_reconstruct_compare.py (řádky 526-533 v něm).
#
#  Pojistky: série B24 nikdy nezapíše do volume_results.csv (B21);
#  ID stromu musí začínat názvem série; chybějící soubor se jen přeskočí.
#
#  Spuštění: nejdřív s DRY_RUN = True (jen vypíše, co by zapsal),
#  teprve potom DRY_RUN = False.
# =====================================================================

import os

import paths
from tree_geom_utils import parse_adqsm_params_file, upsert_result

# =====================  PARAMETRY  ====================================
# Série stromů. Určuje, do kterého souboru se bude zapisovat.
SERIES = "B24"

# Stromy, které se zpracují. ID musí začínat SERIES + "_".
TREE_NAMES = [
    "B24_S01", "B24_S04", "B24_S07", "B24_S08", "B24_S09", "B24_S10",
    "B24_S11", "B24_S12", "B24_S13", "B24_S14", "B24_S16", "B24_S21",
    "B24_S24",
]

# Varianty AdQSM = názvy podsložek stromu. "999" je medián (vyrábí ho
# adqsm_build_median_variant.py); když složka není, přeskočí se.
ADQSM_VARIANTS = ["040", "050", "060", "070", "080", "090", "100", "999"]

# True  = NIC se nezapíše, jen se vypíše tabulka řádků, které by vznikly.
# False = řádky se opravdu zapíšou do souboru s výsledky.
DRY_RUN = True
# =====================================================================


def main():
    # Cílový soubor podle série (neznámá série = chyba z paths.py).
    results_csv = paths.volume_results_csv(SERIES)

    # --- POJISTKA 1: B24 (ani jiná série než B21) nesmí zapisovat do
    # souboru B21. Porovnáváme absolutní cesty.
    b21_csv = os.path.abspath(paths.VOLUME_RESULTS_CSV)
    if SERIES != "B21" and os.path.abspath(results_csv) == b21_csv:
        raise SystemExit(
            "STOP: série %s by zapisovala do %s (soubor série B21). Nic nezapsáno."
            % (SERIES, b21_csv))

    # --- POJISTKA 2: každý strom musí patřit do zvolené série.
    wrong = [t for t in TREE_NAMES if not t.startswith(SERIES + "_")]
    if wrong:
        raise SystemExit(
            "STOP: ID stromů %s nezačínají '%s_' (SERIES = '%s'). Nic nezapsáno."
            % (wrong, SERIES, SERIES))

    print("Série: %s | cílový soubor: %s | DRY_RUN = %s" % (SERIES, results_csv, DRY_RUN))
    print("%-9s %-34s %9s %9s %9s %8s" % ("strom", "metoda", "total_m3", "trunk_m3", "branch_m3", "DBH_m"))

    n_written = 0
    skipped = []   # seznam (strom, varianta, důvod)

    for tree in TREE_NAMES:
        tree_dir = str(paths.tree_data_dir(tree))
        for variant in ADQSM_VARIANTS:
            params_file = os.path.join(tree_dir, variant, "TreesParams.txt")
            if not os.path.exists(params_file):
                skipped.append((tree, variant, "chybí TreesParams.txt"))
                continue

            # Stejné čtení jako v adtree_reconstruct_compare.py (ř. 519).
            # Funkce vrací None, když v souboru není TreeVolume; když chybí
            # TreeDBH, vyhodí chybu - tu zachytíme, ať skript nespadne.
            try:
                adqsm_ref = parse_adqsm_params_file(params_file)
            except Exception as exc:
                skipped.append((tree, variant, "chyba čtení: %s" % exc))
                continue
            if adqsm_ref is None:
                skipped.append((tree, variant, "v TreesParams.txt není TreeVolume"))
                continue

            # Stejný název metody jako v adtree_reconstruct_compare.py
            # (ř. 466 + 527): "AdQSM (TreesParams) (AdQSM 080)".
            method = "AdQSM (TreesParams) (AdQSM %s)" % variant

            print("%-9s %-34s %9.4f %9.4f %9.4f %8.4f"
                  % (tree, method, adqsm_ref["total_vol"], adqsm_ref["trunk_vol"],
                     adqsm_ref["branch_vol"], adqsm_ref["dbh"]))

            if not DRY_RUN:
                # Stejné argumenty jako v adtree_reconstruct_compare.py (ř. 526-533);
                # DBH je v parse_adqsm_params_file() už převedené z cm na m.
                upsert_result(results_csv, tree, method,
                              adqsm_ref["total_vol"], adqsm_ref["trunk_vol"], adqsm_ref["branch_vol"], None,
                              adqsm_ref.get("dbh"), adqsm_ref.get("height"), None,
                              adqsm_ref.get("trunk_len"), adqsm_ref.get("branch_len"),
                              branch_filter="none")
            n_written += 1

    print()
    print("%s: %d řádků, přeskočeno: %d"
          % ("DRY_RUN (nic nezapsáno)" if DRY_RUN else "Zapsáno", n_written, len(skipped)))
    for tree, variant, reason in skipped:
        print("  přeskočeno: %s / %s - %s" % (tree, variant, reason))


if __name__ == "__main__":
    main()
