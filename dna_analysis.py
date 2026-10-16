# =============================================================
#  dna_analysis.py  –  Fonctions d'analyse d'une séquence ADN
#  VERSION CORRIGÉE ET AMÉLIORÉE
# =============================================================
# CORRECTIONS APPORTÉES :
#  1. find_shine_dalgarno : distance étendue à 5-10 nt (était 7-9 → manquait des sites)
#  2. find_all_sd_sites   : NOUVEAU – recherche SD sur les DEUX brins
#  3. find_terminators_both_strands : NOUVEAU – terminateurs sur les deux brins
#  4. Défaut min_polyT abaissé à 3 (était 4 → aucun terminateur trouvé)
#  5. get_gc_content      : NOUVEAU – calcul du contenu GC
# =============================================================

from typing import List, Tuple, Dict, Optional

# ------------------------------------------------------------------
# Table du code génétique (codon → acide aminé)
# ------------------------------------------------------------------
GENETIC_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

STOP_CODONS = {"TAA", "TAG", "TGA"}


# ==================================================================
#  SECTION 1 — Utilitaires de base
# ==================================================================

def clean_sequence(seq: str) -> str:
    """Nettoyer et normaliser la séquence ADN (majuscules, sans espaces)."""
    cleaned = seq.upper().strip()
    for ch in (" ", "\n", "\r", "\t"):
        cleaned = cleaned.replace(ch, "")
    return cleaned


def complement_base(base: str) -> str:
    """Retourner le complément d'une base ADN."""
    return {"A": "T", "T": "A", "G": "C", "C": "G"}.get(base, "N")


def reverse_complement(seq: str) -> str:
    """Calculer le complément inverse (brin antisens 5'→3')."""
    return "".join(complement_base(b) for b in reversed(seq))


def translate(seq: str) -> str:
    """
    Traduire une séquence ADN en protéine.
    La traduction s'arrête au premier codon STOP (non inclus).
    """
    protein = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        aa = GENETIC_CODE.get(codon, "?")
        if aa == "*":
            break
        protein.append(aa)
    return "".join(protein)


def count_mismatches(s1: str, s2: str) -> int:
    """Compter le nombre de positions différentes entre deux séquences."""
    return sum(a != b for a, b in zip(s1, s2))


def get_gc_content(seq: str) -> float:
    """Calculer le pourcentage de GC d'une séquence ADN."""
    if not seq:
        return 0.0
    gc = sum(1 for b in seq if b in "GC")
    return round(gc / len(seq) * 100, 1)


# ==================================================================
#  SECTION 2 — Détection des ORFs (6 cadres de lecture)
# ==================================================================

def find_orfs_in_sequence(seq: str, min_length_aa: int = 50) -> List[Dict]:
    """
    Trouver tous les ORFs complets (ATG→STOP) dans les 3 cadres de lecture
    d'une séquence (brin sens ou antisens).

    Paramètres :
        seq           – séquence ADN 5'→3'
        min_length_aa – longueur minimale en acides aminés (sans le Met initial)

    Retourne une liste de dicts triée par longueur décroissante.
    """
    orfs = []
    for frame in range(3):
        i = frame
        while i < len(seq) - 2:
            codon = seq[i:i + 3]
            if codon == "ATG":
                # Chercher le premier codon STOP en phase, après l'ATG
                found_stop = False
                for j in range(i + 3, len(seq) - 2, 3):
                    stop_codon = seq[j:j + 3]
                    if stop_codon in STOP_CODONS:
                        orf_seq = seq[i:j + 3]
                        aa_len = (j - i) // 3  # nombre d'aa hors Met initial
                        if aa_len >= min_length_aa:
                            orfs.append({
                                "start":      i,
                                "end":        j + 3,
                                "frame":      frame + 1,
                                "length_nt":  len(orf_seq),
                                "length_aa":  aa_len,
                                "sequence":   orf_seq,
                                "protein":    translate(orf_seq),
                                "stop_codon": stop_codon,
                                "gc":         get_gc_content(orf_seq),
                            })
                        i = j + 3
                        found_stop = True
                        break
                if not found_stop:
                    i += 3
            else:
                i += 3

    return sorted(orfs, key=lambda x: x["length_aa"], reverse=True)


def find_all_orfs(seq: str, min_length_aa: int = 50) -> Tuple[List[Dict], List[Dict]]:
    """
    Rechercher les ORFs sur les DEUX brins d'ADN (6 cadres de lecture).

    Retourne :
        sense_orfs     – ORFs du brin sens     (+)
        antisense_orfs – ORFs du brin antisens (−)

    Les coordonnées des ORFs antisens sont converties en positions sur le
    brin sens original (orig_start, orig_end).
    """
    rc_seq = reverse_complement(seq)
    n = len(seq)

    sense_orfs = find_orfs_in_sequence(seq, min_length_aa)
    for orf in sense_orfs:
        orf["strand"] = "+"

    antisense_orfs = find_orfs_in_sequence(rc_seq, min_length_aa)
    for orf in antisense_orfs:
        orf["strand"] = "-"
        # Convertir vers les coordonnées du brin sens original
        orf["orig_start"] = n - orf["end"]
        orf["orig_end"]   = n - orf["start"]

    return sense_orfs, antisense_orfs


# ==================================================================
#  SECTION 3 — Identification de l'ORF codant potentiel
# ==================================================================

def get_best_orf(
    sense_orfs: List[Dict],
    antisense_orfs: List[Dict],
    sd_sites: Optional[List[Dict]] = None
) -> Optional[Dict]:
    """
    Identifier l'ORF le plus probable pour coder une protéine.

    Algorithme de score :
        score = longueur_aa
        + 200  si un site SD valide précède l'ATG
    """
    all_orfs = sense_orfs + antisense_orfs
    if not all_orfs:
        return None

    # Positions ATG associées à un site SD
    sd_atg_positions: set = set()
    if sd_sites:
        sd_atg_positions = {s["atg_pos_orig"] for s in sd_sites}

    best_orf = None
    best_score = -1

    for orf in all_orfs:
        score = orf["length_aa"]
        # Position ATG dans les coordonnées originales
        atg_orig = orf["start"] if orf["strand"] == "+" else orf["orig_start"]
        if atg_orig in sd_atg_positions:
            score += 200
        if score > best_score:
            best_score = score
            best_orf = orf

    return best_orf


# ==================================================================
#  SECTION 4 — Promoteurs bactériens (-35 / -10)
# ==================================================================

def find_promoters(seq: str, max_mismatches: int = 2) -> List[Dict]:
    """
    Rechercher les promoteurs bactériens σ70 avec :
        Box -10 : TATAAT  (consensus)
        Box -35 : TTGACA  (consensus)
        Spacer  : 15–21 nt entre les deux boxes (optimal 17-19)

    Paramètres :
        max_mismatches – mismatches tolérés dans chaque box (0-3)

    Score : 12 = parfait (0 mismatch sur chaque box)
    """
    box10_cons = "TATAAT"
    box35_cons = "TTGACA"
    promoters = []

    for i in range(len(seq) - 6):
        subseq10 = seq[i:i + 6]
        mm10 = count_mismatches(subseq10, box10_cons)
        if mm10 > max_mismatches:
            continue

        for spacer in range(15, 22):          # spacer typique : 17-19 nt
            pos35 = i - spacer - 6
            if pos35 < 0:
                continue
            subseq35 = seq[pos35:pos35 + 6]
            mm35 = count_mismatches(subseq35, box35_cons)
            if mm35 <= max_mismatches:
                promoters.append({
                    "pos_35":  pos35,
                    "seq_35":  subseq35,
                    "mm_35":   mm35,
                    "spacer":  spacer,
                    "pos_10":  i,
                    "seq_10":  subseq10,
                    "mm_10":   mm10,
                    "score":   12 - mm10 - mm35,   # 12 = consensus parfait
                    "optimal_spacer": 17 <= spacer <= 19,
                })

    # Garder le meilleur promoteur par position -10
    seen: Dict = {}
    for p in promoters:
        key = p["pos_10"]
        if key not in seen or p["score"] > seen[key]["score"]:
            seen[key] = p

    return sorted(seen.values(), key=lambda x: x["score"], reverse=True)


# ==================================================================
#  SECTION 5 — Terminateurs de transcription Rho-indépendants
# ==================================================================

def find_palindromes(
    seq: str,
    min_arm: int = 4,
    max_arm: int = 12,
    min_loop: int = 3,
    max_loop: int = 8
) -> List[Dict]:
    """
    Rechercher les séquences palindromiques (structures en épingle à cheveux).

    Un palindrome est de la forme :
        5'—[bras_gauche]—[boucle]—[bras_droit]—3'
    où bras_droit = complément_inverse(bras_gauche).
    """
    palindromes = []
    n = len(seq)

    for i in range(n):
        for arm_len in range(min_arm, max_arm + 1):
            if i + arm_len > n:
                break
            left_arm = seq[i:i + arm_len]
            for loop_len in range(min_loop, max_loop + 1):
                right_start = i + arm_len + loop_len
                if right_start + arm_len > n:
                    break
                right_arm = seq[right_start:right_start + arm_len]
                if left_arm == reverse_complement(right_arm):
                    palindromes.append({
                        "start":      i,
                        "end":        right_start + arm_len,
                        "left_arm":   left_arm,
                        "loop":       seq[i + arm_len:right_start],
                        "right_arm":  right_arm,
                        "arm_length": arm_len,
                        "loop_length":loop_len,
                    })

    # Déduplication : garder le palindrome au bras le plus long par position
    seen: Dict = {}
    for p in palindromes:
        key = p["start"]
        if key not in seen or p["arm_length"] > seen[key]["arm_length"]:
            seen[key] = p

    return list(seen.values())


def find_terminators(seq: str, min_polyT: int = 3) -> List[Dict]:
    """
    Détecter les terminateurs Rho-indépendants (brin sens) :
    palindrome suivi d'une queue poly-T.

    Un terminateur valide : palindrome + ≥ min_polyT 'T' consécutifs en aval.
    """
    palindromes = find_palindromes(seq)
    terminators = []

    for pal in palindromes:
        downstream = seq[pal["end"]:pal["end"] + 15]
        t_count = 0
        for base in downstream:
            if base == "T":
                t_count += 1
            else:
                break
        if t_count >= min_polyT:
            terminators.append({
                **pal,
                "polyT":        downstream[:t_count],
                "polyT_length": t_count,
                "strength":     pal["arm_length"] + t_count,
            })

    return sorted(terminators, key=lambda x: x["strength"], reverse=True)


def find_terminators_both_strands(seq: str, min_polyT: int = 3) -> List[Dict]:
    """
    ✅ NOUVEAU — Rechercher les terminateurs sur les DEUX brins.
    Les terminateurs du brin antisens sont convertis en coordonnées originales.
    """
    rc_seq = reverse_complement(seq)
    n = len(seq)

    sense_terms = find_terminators(seq, min_polyT)
    for t in sense_terms:
        t["strand"] = "+"

    antisense_terms = find_terminators(rc_seq, min_polyT)
    for t in antisense_terms:
        t["strand"] = "-"
        t["orig_start"] = n - t["end"]
        t["orig_end"]   = n - t["start"]

    all_terms = sense_terms + antisense_terms
    return sorted(all_terms, key=lambda x: x["strength"], reverse=True)


# ==================================================================
#  SECTION 6 — Site Shine-Dalgarno (fixation ribosome)
# ==================================================================

def find_shine_dalgarno(seq: str, max_mismatches: int = 2) -> List[Dict]:
    """
    Rechercher les sites Shine-Dalgarno sur une séquence donnée.

    Consensus  : AGGAGG
    Distance   : 5 à 10 nt entre la fin du site SD et le A de l'ATG
                 (✅ CORRIGÉ : était 7-9, étendu à 5-10 pour plus de sensibilité)
    Score max  : 6/6 (consensus parfait)

    Les positions retournées sont relatives à la séquence fournie.
    """
    sd_consensus = "AGGAGG"
    sd_sites = []

    atg_positions = [i for i in range(len(seq) - 2) if seq[i:i + 3] == "ATG"]

    for atg_pos in atg_positions:
        for distance in range(5, 11):          # ✅ 5 à 10 nt (standard élargi)
            sd_start = atg_pos - distance - 6
            if sd_start < 0:
                continue
            subseq = seq[sd_start:sd_start + 6]
            if len(subseq) < 6:
                continue
            mm = count_mismatches(subseq, sd_consensus)
            if mm <= max_mismatches:
                sd_sites.append({
                    "sd_start":  sd_start,
                    "sd_seq":    subseq,
                    "mismatches":mm,
                    "atg_pos":   atg_pos,
                    "distance":  distance,
                    "score":     6 - mm,
                })

    # Déduplication par (position SD, position ATG)
    seen: Dict = {}
    for s in sd_sites:
        key = (s["sd_start"], s["atg_pos"])
        if key not in seen or s["score"] > seen[key]["score"]:
            seen[key] = s

    return sorted(seen.values(), key=lambda x: x["score"], reverse=True)


def find_all_sd_sites(
    seq: str,
    sense_orfs: List[Dict],
    antisense_orfs: List[Dict],
    max_mismatches: int = 2
) -> List[Dict]:
    """
    ✅ NOUVEAU — Rechercher les sites SD sur les DEUX brins.

    - Pour les ORFs sens    : cherche sur le brin sens
    - Pour les ORFs antisens: cherche sur le brin complémentaire inverse

    Retourne tous les sites avec :
        atg_pos_orig – position de l'ATG dans les coordonnées originales
        strand       – '+' ou '-'
    """
    rc_seq = reverse_complement(seq)
    n = len(seq)

    # ATG du brin sens
    sense_atg_pos = {o["start"] for o in sense_orfs}
    # ATG du brin antisens (en coordonnées RC)
    antisense_atg_pos_rc = {o["start"] for o in antisense_orfs}   # dans RC

    # --- Brin sens ---
    sense_sd = []
    for s in find_shine_dalgarno(seq, max_mismatches):
        if s["atg_pos"] in sense_atg_pos:
            s["strand"]       = "+"
            s["atg_pos_orig"] = s["atg_pos"]
            sense_sd.append(s)

    # --- Brin antisens ---
    antisense_sd = []
    for s in find_shine_dalgarno(rc_seq, max_mismatches):
        if s["atg_pos"] in antisense_atg_pos_rc:
            s["strand"]       = "-"
            # Position ATG dans les coords originales
            s["atg_pos_orig"] = n - s["atg_pos"] - 3
            s["sd_start_orig"]= n - s["sd_start"] - 6
            antisense_sd.append(s)

    all_sd = sense_sd + antisense_sd
    return sorted(all_sd, key=lambda x: x["score"], reverse=True)
