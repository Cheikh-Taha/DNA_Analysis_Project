# 🧬 Application d'Analyse d'ADN Bactérien

> Développé par **Taha CHEIKH** — Biologie Moléculaire · Bioinformatique

Application Python interactive permettant l'analyse complète d'une séquence ADN bactérienne : détection des ORFs, promoteurs, terminateurs Rho-indépendants et sites Shine-Dalgarno.

---

## 📋 Table des matières

1. [Description du projet](#description)
2. [Concepts biologiques](#concepts)
3. [Installation](#installation)
4. [Utilisation](#utilisation)
5. [Architecture du code](#architecture)
6. [Fonctions principales](#fonctions)
7. [Résultats obtenus](#resultats)
8. [Paramètres configurables](#parametres)
9. [Bugs corrigés](#bugs)

---

## Description du projet <a name="description"></a>

Ce projet analyse une séquence ADN bactérienne de **1865 nucléotides** (GC = 30.6%, AT = 69.4%) pour en identifier automatiquement les éléments fonctionnels clés :

| Élément                | Rôle biologique                              |
|------------------------|----------------------------------------------|
| ORFs (cadres ouverts)  | Régions potentiellement codantes              |
| Promoteurs –10/–35     | Signal de démarrage de la transcription       |
| Terminateurs Rho-indép.| Signal d'arrêt de la transcription           |
| Sites Shine-Dalgarno   | Signal de démarrage de la traduction         |

L'application est développée avec **Python 3** pour le moteur d'analyse et **Streamlit** pour l'interface web interactive.

---

## Concepts biologiques <a name="concepts"></a>

### 1. ORF — Open Reading Frame

Un **ORF (cadre de lecture ouvert)** est une séquence d'ADN qui :
- Commence par un codon **ATG** (méthionine, codon d'initiation)
- Se termine par un codon **STOP** : `TAA`, `TAG` ou `TGA`
- Est lue dans un des **6 cadres de lecture** possibles

#### Les 6 cadres de lecture

L'ADN est double brin. Pour chaque brin, il y a 3 façons de lire la séquence selon le décalage de départ :

```
Brin sens    (+) : Cadre +1, +2, +3  →  lecture 5'→3' directe
Brin antisens(–) : Cadre –1, –2, –3  →  lecture sur le brin complémentaire inverse
```

Le **brin antisens** est obtenu par la fonction `reverse_complement()` : 
- Complémentarité : A↔T, G↔C
- Inversion : la séquence est lue dans l'autre sens

#### Critère de sélection

Seuls les ORFs d'au moins **50 acides aminés** (≥ 150 nt) sont retenus pour éviter les faux positifs.

#### Choix du meilleur ORF

Un score est attribué selon :
- **Longueur** en acides aminés (principal critère)
- **+200 points** si un site Shine-Dalgarno est associé (dans les 5-10 nt en amont de l'ATG)

---

### 2. Promoteurs bactériens — Box –10 / –35

Un **promoteur** est la région où l'ARN polymérase (avec le facteur σ70) se fixe pour initier la transcription.

Un promoteur σ70 bactérien est composé de deux éléments consensus :

```
5' ... TTGACA ... (17-19 nt) ... TATAAT ... → TSS ... 3'
       Box –35                  Box –10
```

| Élément | Consensus | Rôle                                      |
|---------|-----------|-------------------------------------------|
| Box –35 | `TTGACA`  | Reconnaissance initiale par le facteur σ70 |
| Spacer  | 17–19 nt  | Distance optimale entre les deux boîtes   |
| Box –10 | `TATAAT`  | Ouverture de la double hélice ADN (fusion) |
| TSS     | —         | Site de démarrage de la transcription      |

#### Système de score

Chaque boîte est comparée à son consensus (6 positions). Le score maximum est **12/12** (6 + 6, 0 mismatch). L'algorithme tolère jusqu'à **2 mismatches** par boîte.

---

### 3. Terminateurs Rho-indépendants

Un **terminateur Rho-indépendant** (intrinsèque) est un signal de fin de transcription qui ne nécessite pas la protéine Rho.

#### Structure

```
5' ... [Bras gauche] [Boucle] [Bras droit] ... TTTTT ... 3'
        ←————————— Palindrome (hairpin) —————————→  Queue poly-T
```

La structure se compose de :
1. **Palindrome** : deux séquences inversées complémentaires formant une **épingle à cheveux** (hairpin) dans l'ARNm
2. **Queue poly-T** : au moins 3 thymines consécutives en aval

#### Mécanisme

1. L'ARNm transcrit se replie en épingle → structure thermodynamiquement stable
2. Cette structure bloque physiquement la progression de l'ARN polymérase
3. La queue poly-T déstabilise les liaisons faibles A-U entre l'ARN et le gabarit ADN
4. L'ARN polymérase se détache → fin de la transcription

#### Algorithme de détection

```
1. Chercher des palindromes (bras de 4-12 nt, boucle de 3-8 nt)
2. Vérifier la présence d'au moins 3T en aval (< 10 nt)
3. Calculer une "force" = longueur bras + longueur queue poly-T
4. Appliquer sur les 2 brins (sens + antisens)
```

---

### 4. Site Shine-Dalgarno

Le **site Shine-Dalgarno (SD)** est une séquence courte de l'ARNm qui joue un rôle crucial dans l'initiation de la traduction chez les procaryotes.

#### Principe

La sous-unité 30S du ribosome contient l'**ARNr 16S**, dont l'extrémité 3' est complémentaire du site SD. Cette complémentarité permet au ribosome de se positionner précisément sur le codon ATG d'initiation.

```
ARNm : 5' ... AGGAGG ... (5-10 nt) ... ATG ... 3'
                ↕ complémentarité
ARNr 16S (30S) : 3' ... UCCUCC ... 5'
```

| Paramètre  | Valeur                          |
|------------|---------------------------------|
| Consensus  | `AGGAGG`                        |
| Distance   | 5 à 10 nt en amont de l'ATG    |
| Tolérance  | ≤ 2 mismatches                  |
| Score max  | 6/6 (consensus parfait)         |

---

## Installation <a name="installation"></a>

### Prérequis

- Python 3.8 ou supérieur
- pip

### Installation des dépendances

```bash
pip install streamlit pandas
```

### Cloner le projet

```bash
git clone https://github.com/votre-username/dna-analysis-app.git
cd dna-analysis-app
```

---

## Utilisation <a name="utilisation"></a>

### Lancer l'application

```bash
streamlit run app.py
```

L'application s'ouvre automatiquement dans votre navigateur à l'adresse :
```
http://localhost:8501
```

### Interface

1. **Charger une séquence** : La séquence par défaut (1865 nt) est déjà chargée. Vous pouvez coller votre propre séquence dans le champ texte.
2. **Configurer les paramètres** dans la barre latérale :
   - Longueur minimale des ORFs (en acides aminés)
   - Nombre maximal de mismatches pour les promoteurs
   - Nombre maximal de mismatches pour les sites SD
   - Longueur minimale de la queue poly-T
3. **Cliquer sur "Analyser"**
4. **Consulter les résultats** dans les 5 sections :
   - Section 0 : Statistiques de la séquence
   - Section 1 : Tableau des ORFs (sens + antisens)
   - Section 2 : Meilleur ORF codant + protéine traduite
   - Section 3 : Tableau des promoteurs
   - Section 4 : Tableau des terminateurs
   - Section 5 : Tableau des sites Shine-Dalgarno
   - Récapitulatif final

---

## Architecture du code <a name="architecture"></a>

```
dna_app/
├── dna_analysis.py    # Moteur d'analyse — fonctions biologiques
└── app.py             # Interface Streamlit
```

### Séparation des responsabilités

| Fichier            | Rôle                                        |
|--------------------|---------------------------------------------|
| `dna_analysis.py`  | Logique biologique pure, sans interface     |
| `app.py`           | Affichage, paramètres, mise en forme        |

Cette séparation permet de réutiliser `dna_analysis.py` dans d'autres projets ou scripts sans dépendre de Streamlit.

---

## Fonctions principales <a name="fonctions"></a>

### `dna_analysis.py`

#### Utilitaires

```python
def clean_sequence(seq: str) -> str:
    """Normalise la séquence : majuscules, supprime les espaces/retours à la ligne."""

def reverse_complement(seq: str) -> str:
    """Retourne le brin complémentaire inverse (antisens, 5'→3')."""

def translate(seq: str) -> str:
    """Traduit une séquence ADN en protéine. S'arrête au premier codon STOP."""

def get_gc_content(seq: str) -> float:
    """Retourne le pourcentage de GC de la séquence."""
```

#### Détection des ORFs

```python
def find_orfs_in_sequence(seq, min_length_aa=50):
    """
    Cherche les ORFs dans une séquence (3 cadres).
    Retourne une liste de dicts : {start, end, frame, length_aa, protein, gc, ...}
    """

def find_all_orfs(seq, min_length_aa=50):
    """
    Cherche les ORFs sur les 6 cadres (sens + antisens).
    Retourne : (sense_orfs, antisense_orfs)
    Les coordonnées des ORFs antisens sont converties en coordonnées originales.
    """

def get_best_orf(sense_orfs, antisense_orfs, sd_sites=None):
    """
    Identifie le meilleur ORF codant.
    Score = length_aa + 200 si un site SD est associé dans les 5-10 nt amont.
    """
```

#### Détection des promoteurs

```python
def find_promoters(seq, max_mismatches=2):
    """
    Recherche les promoteurs σ70 sur le brin sens.
    Critères : Box-35 (TTGACA) + spacer 15-21 nt + Box-10 (TATAAT).
    Score de 0 à 12 (6 positions × 2 boîtes).
    """
```

#### Détection des terminateurs

```python
def find_palindromes(seq, min_arm=4, max_arm=12, min_loop=3, max_loop=8):
    """Détecte les structures palindromiques (hairpins) dans la séquence."""

def find_terminators_both_strands(seq, min_polyT=3):
    """
    Détecte les terminateurs Rho-indépendants sur les 2 brins.
    Force = longueur bras + longueur queue poly-T.
    """
```

#### Détection des sites Shine-Dalgarno

```python
def find_shine_dalgarno(seq, max_mismatches=2):
    """
    Recherche le consensus AGGAGG à 5-10 nt d'un ATG.
    Score = nombre de positions correspondant au consensus (sur 6).
    """

def find_all_sd_sites(seq, sense_orfs, antisense_orfs, max_mismatches=2):
    """
    Recherche les sites SD sur les 2 brins.
    Retourne les sites avec leur brin, position ATG originale, et score.
    """
```

---

## Résultats obtenus <a name="resultats"></a>

Analyse de la séquence bactérienne de 1865 nt (GC = 30.6%) :

| Section              | Résultat                                         |
|----------------------|--------------------------------------------------|
| **ORFs (≥50 aa)**    | 2 ORFs : 1 sens (56 aa) + 1 antisens (387 aa)   |
| **Meilleur ORF**     | 387 aa, brin −, cadre 2, pos. 310 (orig.)        |
| **Protéine**         | `MKKQLFILYFNIFLIFLGIG…` (membrane possible)      |
| **Promoteurs**       | 71 détectés, meilleur score : 11/12 (pos. 1759)  |
| **Terminateurs**     | 19 sur les 2 brins, force max = 13 (brin −)      |
| **Sites SD**         | 4 sites, meilleur : AGGGGG score 5/6 (brin −)    |
| **Durée d'exécution**| < 5 secondes                                     |

---

## Paramètres configurables <a name="parametres"></a>

| Paramètre           | Par défaut | Plage    | Description                             |
|---------------------|------------|----------|-----------------------------------------|
| `min_orf`           | 50 aa      | 20–300   | Longueur minimale des ORFs              |
| `max_mm_promoteur`  | 2          | 0–3      | Mismatches max pour les promoteurs      |
| `max_mm_sd`         | 2          | 0–2      | Mismatches max pour les sites SD        |
| `min_polyT`         | 3          | 2–8      | Longueur minimale de la queue poly-T    |

---

## Bugs corrigés <a name="bugs"></a>

| Bug identifié                          | Symptôme                         | Correction appliquée                         |
|----------------------------------------|----------------------------------|----------------------------------------------|
| SD uniquement sur le brin sens         | 0 site SD trouvé                 | Ajout de `find_all_sd_sites()` (2 brins)     |
| Distance SD trop stricte (7-9 nt)      | Site réel à 6 nt manqué          | Plage étendue à 5-10 nt                      |
| Terminateurs sur 1 seul brin           | Terminateurs antisens manqués    | `find_terminators_both_strands()` créée      |
| Seuil poly-T trop élevé (≥4T)         | 0 terminateur détecté            | Seuil abaissé à ≥3T par défaut               |

---

## Auteur

**MAROUANE EL BOUSAIRI**  
Biologie Moléculaire · Bioinformatique

---

*Projet réalisé avec Python 3 · Streamlit · Pandas*
