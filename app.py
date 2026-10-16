# =============================================================
#  app.py  –  Application Streamlit d'analyse d'ADN
#  VERSION CORRIGÉE ET AMÉLIORÉE
#  Lancer avec :  streamlit run app.py
# =============================================================
# CORRECTIONS APPORTÉES :
#  1. SD : utilise find_all_sd_sites (DEUX brins) — était uniquement brin sens
#  2. Terminateurs : utilise find_terminators_both_strands (DEUX brins)
#  3. Valeurs par défaut ajustées (min_orf=50, min_polyT=3, max_mm_sd=2)
#  4. Affichage de la longueur en nt dans les tableaux ORF
#  5. Statistiques GC ajoutées
#  6. Coordonnées antisens correctement affichées
# =============================================================

import streamlit as st
import pandas as pd
from dna_analysis import (
    clean_sequence,
    get_gc_content,
    find_all_orfs,
    get_best_orf,
    find_promoters,
    find_terminators_both_strands,
    find_all_sd_sites,
)

# ── Configuration de la page ──────────────────────────────────────
st.set_page_config(
    page_title="Analyse ADN Bactérien",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── Séquence ADN du sujet ─────────────────────────────────────────
DEFAULT_SEQ = (
    "CATTTCTCTTAAGATTTATTCTATCTTAACACAACAACTTTTAATAAAAGATATGTAGATTACAATTTAAATAGATTGTAAT"
    "ATTTGTAACACTAACATTAATATAGTTGTTATTTTTGTTACATAAACCACTAATAACTCATAATCTTTTAAAACTTATATTT"
    "GAGATAACATCAACTTTACATTACAAGTTATAAAACAAAAGAAGTGGGACACAGAATTCGTCTTGAACACTGTGTCCCACCT"
    "CGTCCCCAAAACTTGCTCTGTCCGTAGAAAAATAAAAAGGGGCCCCCTTTGTTGGGGCCCCGTATTAATTTGCTTCTTTACG"
    "GCGTGACTTAAGTCCTTTTTCAATAAAAATAATGATAATTCCACTTAATGAAACCGCAATAGCCATATATAAAGGAAACTCT"
    "AAATTAACATCGAATAATGCACCAGCTACAAGAGGACCTATAAAGTTACCCATACTGGTAAAAGTTGAATTCAATCCACCTG"
    "CAAAACCTTGCCGTTTGCCTGCTATATTCGAGAAGTAATTGGTTAAAGCTGGTCTAATCATATCAAAACCTATAAAAACAAC"
    "AAAGCTAATAATCATAATCGTCCAATAACCGTTTGCAAGCACTAACATAACGAGAACAATGGCTGAATATAGTAATGACCAT"
    "GCAATAAAATTAAGTTCACTCATATATTTCATAAATTTATCAAAGAAGAATACTTGGAATAATGCGCCAAACACGCCTCCACC"
    "GATAATAGCTATCGAAATATCTTTAGGAGTATAATTTACTTTGTCAGCTGTATATAAAGAAAATAATGTTTCAAAAGCAGATA"
    "AACCAAATGCTAATACAAGTGTTAATATGACTGGAGTAATAAATACTTTCCAATTAATTTTAGTGAATAATTCAGGTTGATAT"
    "TGGTGGAATCCATCTGTAGTTGCTTTTTGAGGATTATGAATTAATAAAACTGACATAATGAATGCAACAACACCTAATGTTCC"
    "AGCAACATAGAAAGGTAATCTATGTGAAATTTCAGCTAAAAAGCCTCCAAATCCAGGTCCTAATATAAAACCTGAATTAATAA"
    "TTGCCGACATGTAACCAAAGTTTTTAGCTTTATCAGCTCCTGGAGAAATATCTGCAATCATACCTGTTACACCAGGCATGACC"
    "ATGCCTGCACTAAAGCCACCTAAAACACGTGAAATGATTAAAATGGTAAAACTTTGACCGGCTGCGAACATAAATTCAGAGACA"
    "GCAAAGAATACTAAACCGATACATATAATTAATTTTTTACCCAATTTATCAGCTAGTGTCCCACCAAATGGTGAAATAATCATT"
    "TGTGATAATGCAAAAGCAGCAACTAGCATTCCTAAGTCACTACCTTTTAATCCTAAATCCTTCAAATATACAGGAAGTACAGGA"
    "ATAACTAATCCAATCCCTAAAAATATAAGAAATATATTAAAATAAAGAATGAATAATTGTTTTTTCATAGGAACCCCCCTTTTT"
    "TTCTAATATAGTAAGCATACCAATAATTTTGAAATTTATCGCAAATTTATTTTTAATGTCCCACATAAGAACTTTACTTTACTA"
    "AAAGTTTCTACATTGTATAACAATTTAACGTGACATGACAATCATTTTGAAAAACTTTAATTTTTAAATTTCACTAAACTCTAAT"
    "TTTTTAGTTAAAAATTTAACATATAATTAAACGAACCTGTATCTTCATAAACTTTCGACAAAACCTTTATCCAATTGACTAAAGT"
    "CTAATGAACATACGGTATAATTTTTCTTCGCTTTCACATCTGTATAAACAGAAAACACCTCAAATTGTTTGAAATGCCTCACTGT"
    "AAAAGTTATCATTAACTCTCGATTATCAAACTCTAT"
)

# ==================================================================
#  SIDEBAR — Paramètres de l'analyse
# ==================================================================
with st.sidebar:
    st.title("⚙️ Paramètres")
    st.markdown("---")

    min_orf = st.slider(
        "Longueur minimale ORF (aa)", 20, 300, 50, step=10,
        help="Taille minimale d'un ORF en acides aminés. Réduire pour voir plus d'ORFs."
    )
    max_mm_promoter = st.slider(
        "Mismatches max – Promoteur", 0, 3, 2,
        help="Différences tolérées par rapport au consensus TATAAT / TTGACA"
    )
    max_mm_sd = st.slider(
        "Mismatches max – Shine-Dalgarno", 0, 2, 2,
        help="Différences tolérées par rapport au consensus AGGAGG"
    )
    min_polyT = st.slider(
        "polyT min – Terminateur", 2, 8, 3,
        help="Nombre minimum de T consécutifs après le palindrome"
    )

    st.markdown("---")
    st.markdown("**ℹ️ Corrections v2**")
    st.caption("✅ SD sur les 2 brins")
    st.caption("✅ Terminateurs sur les 2 brins")
    st.caption("✅ Distance SD étendue 5-10 nt")
    st.caption("✅ polyT défaut abaissé à 3")

# ==================================================================
#  PAGE PRINCIPALE
# ==================================================================
st.title("🧬 Application d'Analyse d'ADN Bactérien")
st.markdown(
    "Détection des **ORFs** (6 cadres), **promoteurs** σ70, "
    "**terminateurs** Rho-indépendants et sites **Shine-Dalgarno** — sur les **deux brins**."
)
st.markdown("---")

# Zone de saisie
seq_input = st.text_area(
    "📋 Entrez ou collez votre séquence ADN :",
    value=DEFAULT_SEQ,
    height=130,
    help="Séquence ADN 5'→3', lettres A/T/G/C."
)

col_btn, col_info = st.columns([2, 5])
with col_btn:
    run = st.button("🔬 Lancer l'analyse", type="primary", use_container_width=True)

# ==================================================================
#  BLOC PRINCIPAL — exécuté au clic du bouton
# ==================================================================
if run:
    seq = clean_sequence(seq_input)

    if not seq or len(seq) < 30:
        st.error("⚠️ Séquence trop courte ou invalide. Minimum 30 nucléotides.")
        st.stop()

    with col_info:
        gc = get_gc_content(seq)
        at = round(100 - gc, 1)
        st.info(f"Séquence : **{len(seq)} nt**  |  GC = **{gc}%**  |  AT = **{at}%**")

    progress = st.progress(0, text="Initialisation…")

    # ================================================================
    #  SECTION 0 — Statistiques de la séquence
    # ================================================================
    st.markdown("---")
    st.header("📊  Statistiques de la Séquence")

    col_a, col_b, col_c, col_d, col_e = st.columns(5)
    col_a.metric("Longueur", f"{len(seq)} nt")
    col_b.metric("Contenu GC", f"{gc}%")
    col_c.metric("Contenu AT", f"{at}%")
    col_d.metric("A", f"{seq.count('A')}")
    col_e.metric("T", f"{seq.count('T')}")

    col_f, col_g, col_h, col_i, col_j = st.columns(5)
    col_f.metric("G", f"{seq.count('G')}")
    col_g.metric("C", f"{seq.count('C')}")
    col_h.metric("Codons ATG (sens)", f"{seq.count('ATG')}")
    from dna_analysis import reverse_complement
    rc_seq = reverse_complement(seq)
    col_i.metric("Codons ATG (antisens)", f"{rc_seq.count('ATG')}")
    col_j.metric("Brins analysés", "2 (+/−)")

    # ================================================================
    #  SECTION 1 — ORFs sur les 6 cadres de lecture
    # ================================================================
    progress.progress(10, text="Recherche des ORFs sur les 6 cadres…")
    sense_orfs, antisense_orfs = find_all_orfs(seq, min_orf)
    total_orfs = len(sense_orfs) + len(antisense_orfs)

    st.markdown("---")
    st.header("1️⃣  Détection des ORFs (6 cadres de lecture)")
    st.caption(
        f"Brin sens : cadres +1, +2, +3  |  "
        f"Brin antisens : cadres −1, −2, −3  |  "
        f"Seuil : ≥ {min_orf} acides aminés"
    )

    m1, m2, m3 = st.columns(3)
    m1.metric("ORFs brin sens (+)", len(sense_orfs))
    m2.metric("ORFs brin antisens (−)", len(antisense_orfs))
    m3.metric("Total", total_orfs)

    tab_s, tab_a = st.tabs(
        [f"➕ Brin Sens ({len(sense_orfs)} ORF)", f"➖ Brin Antisens ({len(antisense_orfs)} ORF)"]
    )

    def make_orf_table(orfs, strand):
        rows = []
        for o in orfs:
            start = o["start"] if strand == "+" else o["orig_start"]
            end   = o["end"]   if strand == "+" else o["orig_end"]
            rows.append({
                "Cadre":        f"{strand}{o['frame']}",
                "Début (nt)":   start,
                "Fin (nt)":     end,
                "Long. (nt)":   o["length_nt"],
                "Long. (aa)":   o["length_aa"],
                "Codon STOP":   o["stop_codon"],
                "GC %":         o["gc"],
                "Protéine":     o["protein"][:25] + ("…" if len(o["protein"]) > 25 else ""),
            })
        return pd.DataFrame(rows)

    with tab_s:
        if not sense_orfs:
            st.warning(f"Aucun ORF ≥ {min_orf} aa sur le brin sens. Réduisez le seuil dans la sidebar.")
        else:
            st.dataframe(make_orf_table(sense_orfs, "+"), use_container_width=True, hide_index=True)
            for i, orf in enumerate(sense_orfs[:3]):
                with st.expander(f"🔍 ORF sens {i+1} — cadre +{orf['frame']} — {orf['length_aa']} aa — pos {orf['start']}"):
                    c1, c2 = st.columns(2)
                    with c1:
                        st.markdown("**Séquence nucléotidique (60 premiers nt) :**")
                        st.code(orf["sequence"][:60] + "…")
                    with c2:
                        st.markdown("**Séquence protéique (30 premiers aa) :**")
                        st.code(orf["protein"][:30] + "…")

    with tab_a:
        if not antisense_orfs:
            st.warning(f"Aucun ORF ≥ {min_orf} aa sur le brin antisens. Réduisez le seuil.")
        else:
            st.dataframe(make_orf_table(antisense_orfs, "−"), use_container_width=True, hide_index=True)
            for i, orf in enumerate(antisense_orfs[:3]):
                with st.expander(
                    f"🔍 ORF antisens {i+1} — cadre −{orf['frame']} — "
                    f"{orf['length_aa']} aa — pos orig. {orf['orig_start']}"
                ):
                    c1, c2 = st.columns(2)
                    with c1:
                        st.markdown("**Séquence (brin antisens, 60 premiers nt) :**")
                        st.code(orf["sequence"][:60] + "…")
                    with c2:
                        st.markdown("**Séquence protéique (30 premiers aa) :**")
                        st.code(orf["protein"][:30] + "…")

    # ================================================================
    #  SECTION 2 — ORF codant potentiel
    # ================================================================
    progress.progress(30, text="Identification de l'ORF codant…")

    # SD calculé ici car nécessaire pour le score de l'ORF codant
    sd_sites_all = find_all_sd_sites(seq, sense_orfs, antisense_orfs, max_mm_sd)
    best_orf = get_best_orf(sense_orfs, antisense_orfs, sd_sites_all)

    st.markdown("---")
    st.header("2️⃣  ORF Potentiellement Codant pour une Protéine")
    st.caption(
        "Critère : ORF le plus long (tous cadres et brins confondus) + bonus +200 pts si site SD associé"
    )

    if best_orf:
        strand_b = best_orf["strand"]
        start_b  = best_orf["start"] if strand_b == "+" else best_orf["orig_start"]
        end_b    = best_orf["end"]   if strand_b == "+" else best_orf["orig_end"]
        frame_b  = best_orf["frame"]

        st.success(
            f"✅ Meilleur ORF identifié — Brin **{strand_b}**, "
            f"cadre **{strand_b}{frame_b}**, position **{start_b}–{end_b}**, "
            f"**{best_orf['length_aa']} acides aminés**"
        )

        b1, b2, b3, b4, b5 = st.columns(5)
        b1.metric("Position début", start_b)
        b2.metric("Position fin", end_b)
        b3.metric("Longueur (aa)", best_orf["length_aa"])
        b4.metric("Longueur (nt)", best_orf["length_nt"])
        b5.metric("Brin / Cadre", f"{strand_b}{frame_b}")

        with st.expander("📜 Séquence protéique complète (acides aminés)"):
            protein = best_orf["protein"]
            # Affichage par blocs de 10 aa
            blocks = [protein[i:i+10] for i in range(0, len(protein), 10)]
            lines  = ["  ".join(blocks[i:i+6]) for i in range(0, len(blocks), 6)]
            st.code("\n".join(lines))
            st.caption(f"Longueur totale : {len(protein)} acides aminés | GC de l'ORF : {best_orf['gc']}%")

        with st.expander("🧬 Séquence nucléotidique de l'ORF"):
            nuc = best_orf["sequence"]
            nuc_lines = [nuc[i:i+60] for i in range(0, len(nuc), 60)]
            for idx, line in enumerate(nuc_lines):
                st.code(f"{idx*60+1:6d}  {line}")
            st.caption(f"Longueur : {len(nuc)} nt | Commence par ATG, finit par {best_orf['stop_codon']}")
    else:
        st.warning("Aucun ORF codant identifié avec les paramètres actuels. Réduisez la longueur minimale.")

    # ================================================================
    #  SECTION 3 — Promoteurs bactériens
    # ================================================================
    progress.progress(50, text="Recherche des promoteurs…")
    promoters = find_promoters(seq, max_mm_promoter)

    st.markdown("---")
    st.header("3️⃣  Promoteurs Bactériens (Box -10 / Box -35)")

    st.info(
        "**Box -10 consensus :** `TATAAT`  |  "
        "**Box -35 consensus :** `TTGACA`  |  "
        "**Spacer optimal :** 17–19 nt  |  "
        "**Score max :** 12/12"
    )

    if promoters:
        st.success(f"✅ {len(promoters)} promoteur(s) détecté(s) — Top 10 affichés")

        promo_rows = []
        for p in promoters[:20]:
            promo_rows.append({
                "Pos. -35":    p["pos_35"],
                "Séq. -35":   p["seq_35"],
                "MM -35":     p["mm_35"],
                "Spacer (nt)":p["spacer"],
                "Optimal?":   "✅" if p["optimal_spacer"] else "⚠️",
                "Pos. -10":   p["pos_10"],
                "Séq. -10":   p["seq_10"],
                "MM -10":     p["mm_10"],
                "Score (/12)":p["score"],
            })
        st.dataframe(pd.DataFrame(promo_rows), use_container_width=True, hide_index=True)

        for i, p in enumerate(promoters[:5]):
            spacer_note = "✅ spacer optimal" if p["optimal_spacer"] else f"⚠️ spacer non-optimal ({p['spacer']} nt)"
            with st.expander(
                f"🔍 Promoteur {i+1} — Score {p['score']}/12 — "
                f"Box -35 @ {p['pos_35']} — Box -10 @ {p['pos_10']} — {spacer_note}"
            ):
                ctx_s = max(0, p["pos_35"] - 5)
                ctx_e = min(len(seq), p["pos_10"] + 15)
                context = seq[ctx_s:ctx_e]
                # Marquer les boxes dans le contexte
                st.markdown(f"""
| Élément    | Position   | Séquence     | Mismatches |
|------------|-----------|--------------|------------|
| **Box -35**| {p['pos_35']} | `{p['seq_35']}` | {p['mm_35']} |
| Spacer     | —         | —            | {p['spacer']} nt |
| **Box -10**| {p['pos_10']} | `{p['seq_10']}` | {p['mm_10']} |
                """)
                st.markdown("**Contexte génomique :**")
                st.code(context)
    else:
        st.warning("Aucun promoteur détecté. Augmentez le nombre de mismatches tolérés (sidebar).")

    # ================================================================
    #  SECTION 4 — Terminateurs Rho-indépendants (deux brins)
    # ================================================================
    progress.progress(70, text="Recherche des terminateurs (2 brins)…")
    terminators = find_terminators_both_strands(seq, min_polyT)

    st.markdown("---")
    st.header("4️⃣  Terminateurs de Transcription Rho-Indépendants")

    st.info(
        "Structure : **palindrome** (épingle à cheveux) suivi d'une **queue poly-T**  |  "
        "Recherche effectuée sur les **deux brins** ✅"
    )

    if terminators:
        st.success(f"✅ {len(terminators)} terminateur(s) détecté(s) sur les deux brins")

        term_rows = []
        for t in terminators[:20]:
            pos  = t.get("orig_start", t["start"]) if t["strand"] == "-" else t["start"]
            pos2 = t.get("orig_end",   t["end"])   if t["strand"] == "-" else t["end"]
            term_rows.append({
                "Brin":          t["strand"],
                "Début (orig.)": pos,
                "Fin (orig.)":   pos2,
                "Bras (nt)":     t["arm_length"],
                "Boucle":        t["loop"],
                "Long. boucle":  t["loop_length"],
                "Queue polyT":   t["polyT"],
                "Long. polyT":   t["polyT_length"],
                "Force":         t["strength"],
            })
        st.dataframe(pd.DataFrame(term_rows), use_container_width=True, hide_index=True)

        for i, t in enumerate(terminators[:5]):
            pos  = t.get("orig_start", t["start"]) if t["strand"] == "-" else t["start"]
            with st.expander(
                f"🔍 Terminateur {i+1} — Brin {t['strand']} — "
                f"pos {pos} — Force {t['strength']}"
            ):
                st.markdown(f"""
**Structure en épingle à cheveux :**

`5'—[ {t['left_arm']} ]—({t['loop']})—[ {t['right_arm']} ]—{t['polyT']}…—3'`

| Partie          | Séquence           | Longueur    |
|-----------------|--------------------|-------------|
| Bras gauche     | `{t['left_arm']}`  | {t['arm_length']} nt |
| Boucle          | `{t['loop']}`      | {t['loop_length']} nt |
| Bras droit      | `{t['right_arm']}` | {t['arm_length']} nt |
| Queue poly-T    | `{t['polyT']}`     | {t['polyT_length']} nt |
| **Force totale**|                    | **{t['strength']}** |
                """)
    else:
        st.warning(
            "Aucun terminateur détecté. Réduisez le polyT minimum (sidebar) à 2 ou 3."
        )

    # ================================================================
    #  SECTION 5 — Sites Shine-Dalgarno (deux brins)
    # ================================================================
    progress.progress(90, text="Recherche des sites Shine-Dalgarno (2 brins)…")

    st.markdown("---")
    st.header("5️⃣  Sites Shine-Dalgarno (Fixation du Ribosome)")

    st.info(
        "Consensus : `AGGAGG`  |  "
        "Distance : **5 à 10 nt** avant l'ATG ✅  |  "
        "Recherche sur les **deux brins** ✅"
    )

    if sd_sites_all:
        st.success(f"✅ {sd_sites_all[0]['strand']} premier site — {len(sd_sites_all)} site(s) Shine-Dalgarno détecté(s)")

        sd_rows = []
        for s in sd_sites_all:
            sd_rows.append({
                "Brin":          s["strand"],
                "Pos. SD (orig)":s.get("sd_start_orig", s["sd_start"]) if s["strand"] == "-" else s["sd_start"],
                "Séquence SD":   s["sd_seq"],
                "Mismatches":    s["mismatches"],
                "Pos. ATG (orig)":s["atg_pos_orig"],
                "Distance (nt)": s["distance"],
                "Score (/6)":    s["score"],
            })
        st.dataframe(pd.DataFrame(sd_rows), use_container_width=True, hide_index=True)

        for i, s in enumerate(sd_sites_all[:5]):
            with st.expander(
                f"🔍 Site SD {i+1} — Brin {s['strand']} — "
                f"Score {s['score']}/6 — ATG @ pos orig. {s['atg_pos_orig']}"
            ):
                st.markdown(f"""
| Élément       | Position (orig.)                                         | Séquence       |
|---------------|----------------------------------------------------------|----------------|
| Site SD       | {s.get('sd_start_orig', s['sd_start']) if s['strand']=='-' else s['sd_start']} | `{s['sd_seq']}` |
| Spacer        | —                                                        | {s['distance']} nt |
| Codon ATG     | {s['atg_pos_orig']}                                      | `ATG`          |
                """)
                # Contexte sur le bon brin
                if s["strand"] == "+":
                    ctx = seq[max(0, s["sd_start"] - 3): s["atg_pos"] + 6]
                    st.markdown("**Contexte brin sens :**")
                else:
                    from dna_analysis import reverse_complement
                    rc = reverse_complement(seq)
                    ctx = rc[max(0, s["sd_start"] - 3): s["atg_pos"] + 6]
                    st.markdown("**Contexte brin antisens (RC) :**")
                st.code(ctx)
    else:
        st.warning(
            "Aucun site SD détecté. Augmentez les mismatches max SD (sidebar) à 2."
        )

    # ================================================================
    #  FIN — Résumé
    # ================================================================
    progress.progress(100, text="Analyse terminée ✅")
    st.markdown("---")

    st.success("🎉 Analyse complète terminée !")

    # Tableau récapitulatif
    st.subheader("📋 Résumé")
    recap = pd.DataFrame([
        {"Section": "1. ORFs",         "Résultat": f"{total_orfs} ORF(s) trouvé(s) sur les 6 cadres"},
        {"Section": "2. ORF codant",   "Résultat": f"{best_orf['length_aa']} aa — brin {best_orf['strand']} cadre {best_orf['frame']}" if best_orf else "Non trouvé"},
        {"Section": "3. Promoteurs",   "Résultat": f"{len(promoters)} promoteur(s) | Meilleur score {promoters[0]['score']}/12" if promoters else "Aucun"},
        {"Section": "4. Terminateurs", "Résultat": f"{len(terminators)} terminateur(s) sur les 2 brins" if terminators else "Aucun"},
        {"Section": "5. Shine-Dalgarno","Résultat": f"{len(sd_sites_all)} site(s) SD sur les 2 brins" if sd_sites_all else "Aucun"},
    ])
    st.dataframe(recap, use_container_width=True, hide_index=True)
