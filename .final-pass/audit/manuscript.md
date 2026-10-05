# Audit: docs/manuscript.tex (branch figures-v2), read-only

Build used for verification: pdflatex, then bibtex, then pdflatex twice, into /tmp/msaudit. Result: 25 pages, 0 overfull and
0 underfull boxes, 0 undefined references or citations. There are two warnings: "Float too large for page by 16.62926pt on input line 458"
and "caption: \setcaptiontype outside box or environment on input line 897".
Constraint on any fix: build.py `check_numbering` (build.py L41-48) reads labels named `fig:<n>` per `\subsection` and requires
exactly 12 occurrences of `\caption{`. Keep the names fig:1..fig:12. Never title a table with `\caption{` (use `\caption*{`,
`\captionof{` or a run-in head).
Line numbers refer to docs/manuscript.tex at HEAD 4f2626a.

## 1. Preamble and macros

Status (no action): all 11 packages are used. amsmath/amssymb (L4) provide pmatrix, aligned, \operatorname, \tfrac, \mathbb, \mathfrak,
\nmid and \varnothing. geometry is L5. fontenc is L6. graphicx and \graphicspath are L7-8, and both paths resolve. caption (L9) serves \caption*, \captionsetup L897 and
the options. placeins[section] (L10) is used: the hook is still active on \section and \section* after titlesec, checked via \meaning\section.
etoolbox (L11) serves \AtBeginEnvironment at L33. array (L12) serves `>{}` at L1072 and L1096. longtable (L13) serves L1096. titlesec is L14, used at L22-24. hyperref
(L15) serves \ref and the links. All 13 macros are used and there are no dead macros: \R, \conf, \hilb, \hamb, \spin, \act, \stat, \X, \Pmat, \Rmat, \Uop,
\Pproj and \jdens, plus the local \closingtable (L904, used 4 times). The following typed forms do not occur outside L41-53: \mathbb{R}, \mathbb R, \mathcal{C},
\mathcal{H}, \mathbf{X}, \mathbf{R}, \mathbf{P}, \mathsf{P}, \chi_{\mathrm{stat}} and j_\varphi.

- L10: the comment "floats stay within their section" is inaccurate. The barrier fires at \section, which here means parts I-IV and the \section* headings. It does not fire at the numbered sections 1-12, which are \subsection (see L17). Fix: "% floats stay within their part (I-IV)".
- L11: etoolbox is loaded for one call (L33). Fix: replace with `\AddToHook{env/itemize/begin}{\setlength{\itemsep}{3pt}\setlength{\parskip}{0pt}}` and drop etoolbox.
- L13: longtable is used only for the Notation table, which fits on one page (p. 24) and has no `\endhead`. Fix: set Notation as a tabular through the same macro as Table D and drop longtable. Otherwise keep it and add `\endhead` so a page break repeats nothing wrong.
- L15: hyperref leaves the PDF Title and Author empty (pdfinfo). The \section* headings (L895, L1094, References) get no bookmarks: manuscript.out lists only I-IV and 1-12. Fix: `\usepackage[hidelinks,pdftitle={Quantum State Spaces of Nonrigid Molecules and their Complexes},pdfauthor={Isaac AshLind}]{hyperref}`, and `\pdfbookmark[1]{<title>}{<key>}` before each \section*.
- L6: T1 fontenc without lmodern depends on cm-super being installed. It is installed here, and SFRM Type 1 fonts are embedded. Fix: add `\usepackage{lmodern}` after L6 for portable Type 1 fonts.
- L46: \act is bypassed by a bare `\cdot` for the group action at L249 (twice), L291 and L295. The output is identical because \cdot is already \mathbin, so this is source consistency only. Fix: `\act` at all four. Keep `p\cdot q` at L838 and L1005, which is a dot product.
- L51: \Uop is bypassed at L719 by `\widetilde{\mathsf U}_s`. Fix: add `\newcommand{\Uopt}[1]{\widetilde{\mathsf U}_{#1}}` next to L51 and use it at L719.
- L49-52: the macro interface is mixed. \Pmat{s}, \Uop{s} and \Pproj{\Gamma} take the subscript as an argument, while \Rmat and \X take a bare subscript (`\Rmat_h`, `\X_0`). Fix: optional. Either give \Rmat an optional subscript argument, or state the two forms in the L40 comment.
- L904: \closingtable is defined with \newcommand in the body, inside the L903 group. Fix: move it to the preamble as the single table macro (see section 2, C1).
- No macros exist for \mathcal F (29 uses), \mathcal U (11) or \mathbb C (14), although \conf and \R exist for the parallel symbols. Fix: optional. Add `\newcommand{\feas}{\mathcal F}`, `\newcommand{\cell}{\mathcal U}` and `\newcommand{\C}{\mathbb C}`, or leave them typed (they are typed uniformly).

## 2. Tables

Inventory, with natural widths measured at the settings in force, against \linewidth = 469.76pt:
- Section 8 weights table, L586-596: `center` and `\small`, `tabular{lcc}`, three \hline, default \arraystretch and \tabcolsep, numbers in text mode, no title, introduced by the colon at L584-585. Width 210.4pt.
- Table A, L1015-1026: \closingtable, which applies \footnotesize, \arraystretch 1.05, \tabcolsep 5pt, a minipage, a centred bold run-in title above, and `\\[3pt]`. Spec `@{}lllrrrrrrrrr@{}`. Width 426.0pt.
- Table B, L1028-1039: `@{}lcccclccc@{}`. Width 362.9pt.
- Table C, L1041-1052: `@{}clccc@{}`. Width 181.5pt.
- Table D, L1071-1086: `@{}p{6.6cm}p{9.2cm}@{}` (raggedright). Width 459.6pt.
- Notation, L1095-1127: `{\small}` longtable `@{}p{5.9cm}p{8.4cm}@{}`, default \tabcolsep 6pt, \hline top and bottom only, no header row. Width 418.9pt, centred by longtable.
All six fit their width, and the log has no overfull box. Rules are uniform: one \hline above the header, one below it and one at the bottom, with no booktabs and no double rules.
Species cells are all `\mathrm`, sums use `\oplus` throughout, and no table cell contains \times or \cong.

Proposed canonical format:
- C1: one preamble macro for every table, promoted from \closingtable: `\refstepcounter{table}\label{#1}`, with `\renewcommand{\thetable}{\Alph{table}}`. The Section 8 table then becomes Table A and the methane tables B-E, and every citation becomes `Table~\ref{...}`. Each table gets a run-in title `\textbf{Table~\thetable.}` above it, in the caption font (\small), justified (only the tabular centred), with fixed skips in place of `\\[3pt]` and `\medskip`.
- C2: body in `\small` with `\arraystretch` 1.05 and `\tabcolsep` 5pt, and `@{}` at both outer edges. Table A measures 453.4pt at \small, so it fits.
- C3: rules stay as they are: \hline above the header, below the header and at the bottom.
- C4: columns. Label, species and group columns are `l`. Integer columns are `r`, with every number in math mode. The two text tables (D, Notation) use `p{\dimexpr.4\linewidth-\tabcolsep}` and `p{\dimexpr.6\linewidth-\tabcolsep}`, which sum to \linewidth.
- C5: headers are lowercase words or symbols, with one word per quantity: `$\Gamma$` for the species column, and `even` and `odd` for the parities.
- C6: caption labels use a period throughout, set with `labelsep=period` in L9 (see the caption-label item below).

Deviations from that format:
- L586-587: the weights table uses `center` and `\small` with spec `lcc`. It has no `@{}` and default \arraystretch and \tabcolsep. Fix: C1 and C2, spec `@{}lrr@{}`.
- L591-593: the weights table has numbers in text mode (`36`, `12`, `24`), while Tables A-C use math mode (`$5$`). Fix: `$36$`, `$12$`, `$24$` and so on.
- L587, L1029, L1042: integer columns are `c`, while Table A (L1016) uses `r`. Fix: weights `@{}lrr@{}`, Table B `@{}lrrrllrrl@{}`, Table C `@{}rlrrr@{}`.
- L1029: in Table B, column 5 ("on $T$", species) is `c` and column 9 ($T/\ker\Gamma$, group names) is `c`, while column 6 (pairs) and Table C column 2 (species) are `l`. Fix: `l` for species and group columns (see the previous item).
- L589: the weights header reads "spatial species & even parity & odd parity", while Table B (L1031) and Table C (L1044) use `$\Gamma$`/`even`/`odd`. Fix: `$\Gamma$ & even & odd`.
- L586: the weights table has no title and cannot be referenced. It is introduced only by the colon at L585. Fix: C1 (run-in title in the author's words, plus \label).
- L1072 and L1096: the two text tables have different column splits. Table D is 6.6+9.2 cm (459.6pt) and Notation is 5.9+8.4 cm (418.9pt, indented about 25pt on each side). Neither equals \linewidth. Fix: C4 for both.
- L1095-1096: Notation uses `\small` with default \tabcolsep 6pt and \arraystretch 1.0, while \closingtable uses 5pt and 1.05. Fix: C2.
- L904: Tables A-D use `\footnotesize` for body and title, while the weights table, the Notation table and the figure captions (L9 `font=small`) are small. Fix: C2 (`\small`).
- L904 (rendered p. 22): `\centering` sits inside the minipage, so Table B's three-line title (L1028) is centred line by line, while multi-line figure captions are justified. Fix: apply `\centering` to the tabular only.
- L9 against L1015, L1028, L1041 and L1071: figure captions render "Figure 7:" (caption's default colon), while table titles render "Table A.". Fix: `\usepackage[font=small,labelfont=bf,labelsep=period]{caption}`.
- L904: the title gap `\\[3pt]` and the `\par\medskip` around each table are manual skips, and the weights table gets `center`'s \topsep instead. Fix: C1 (one macro, one pair of skips).
- L1020: in Table A's `$\Rmat_h$` column, row E holds math `$1$` while every other row holds words. Fix: `identity`.
- L1015-1052: Tables A-C sit after step 12, while their first citations are L925 (A, p. 20), L958 (B) and L1013 (C). Table D (L1071) follows its citation at L1060. Fix: none needed if this is the intended layout (L892 "each table is kept whole"). Otherwise place each table after the paragraph that first cites it.
- L1028 and L1074: Table B's title says "(their Table~II, $T$ block)" and Table D's header says "Albert et al.", neither with a citation. Albert et al. are first named in the closing section at L1054-1055, after Table B. Fix: `Albert et al.~\cite{AlbertEtAlOrientation}` in the L1074 header, and the same in place of "their" in L1028.
- L1096-1097: Notation has no header row, unlike Table D, the other two-column text table (header at L1074). Fix: optional. Add `symbol & meaning\\ \hline` after L1097.
- L1110, L1118, L1119, L1124, L1125: the Notation separator convention is broken. The convention is ";" between symbols with separate glosses and "," between symbols that share one gloss, as in L1098, L1107, L1121 and L1123. L1110 (`$\mathbf A_s(a)$; $\mathbf M_s(a)$`) and L1118 (`$\stat$; $\chi_\pm$`) use ";" with a shared gloss. L1119, L1124 and L1125 use "," with separate glosses. Fix: "," at L1110 and L1118. At L1119, L1124 and L1125 use ";" between the symbols and between the glosses, for example `$\Gamma$; $\chi_\Gamma$; $d_\Gamma$; $w_\Gamma$ & species; character; dimension; projection share`.

## 3. Notation consistency

Status (verified uniform, no action):
- Species are `\mathrm A_1`, `\mathrm A_2`, `\mathrm E`, `\mathrm T_{1,2}`, `\mathrm A`, `\mathrm T`, `{}^1\mathrm E` and `{}^2\mathrm E` everywhere: text, captions, Tables A-C, Notation, 124 occurrences (one split across L948-949), with no `A_1`, `\text{A}` or `\textrm` variants.
- Group names: `G_6` (15 occurrences), `G_{12}` (10), `SO(3)` (24, italic throughout) and `T_d(\mathrm M)` (2).
- `\mathbb C` is always unbraced in the body (14 occurrences). `[H,B]` and `[H,S]` never have a space (L425, L448, L455-456, L1102).
- `\omega` is used only for angles (L701, L744-752, L853, L983-984, L1010, L1083, L1109, L1122). `\eta` is used only for packets (L612-621, L961-964, L1079, L1121). `\iota` is used only for the umbrella coordinate (L476-477, L504, L730-760, L773, L831, L1112). `\Delta` is used only for packet width (L630-647, L1121); "position variance $\Delta^2$" (L630) agrees with "density standard deviation" (L1121).
- `q` is always the normal coordinates with `n_q` (L665, L831, L839, L843); for methane n_q = 9 (L936, L974, L983, L1005).
- `\ref` and quote marks are covered in section 4.

Inconsistencies:
- L995: `O/T\cong\mathbb Z_2`. Elsewhere cyclic groups are written `C_n`: `C_2` at L112 and `C_2^3` at L321, `C_3` at L1000, L1035 and L1083. `\mathbb Z` means the integers at L842 and L857. Fix: `O/T\cong C_2`.
- L942: `\mathbb C[G/H]=\mathrm A_1`, while the same kind of statement uses `\cong` at L541, L556 and L946. Fix: `\cong`.
- L1124: Notation glosses `\ell` only as the "water ... bond scale (figure 1)". In the text, `\ell` appears only at L912, as methane's C-H length (1.087 Å). Figure 1's caption (L156-159) does not mention `\ell`, although the plate labels "ℓ±δ". Fix: gloss "bond scale (water, Figure~\ref{fig:1}; C--H, closing example)".
- L1106: Notation has `$\gamma_s$` with "$\X\to s\X$". The text introduces an unsubscripted `\gamma` (L387) and writes the general path as `\gamma_g` from `\X` to `g\X` at L793. The figure loops `\gamma_t`, `\gamma_b` (L819) and `\gamma_g` (L997) are specific cases. Fix: L793 `$\gamma_s$ from $\X$ to $s\X$` ($s\in G$), matching L1106, and keep the figure loops.
- L1098: the Notation cell reads "centred configuration space", while L87 reads "centered configuration space". The comment at L42 also says "centred". Fix: one spelling in the table and the comment, matching the text.
- The letter T is overloaded. It denotes the tetrahedral group (L930 onward, Tables B and D), the coordinate action `T_s`, `T_h`, `T_g` (L684, L962, L978, L986-989), the tangent space `T_{\X_0}\conf` (L948), the transpose `^{\mathsf T}` (L76, L86) and the species `\mathrm T`. Group and action meet in one passage at L986-989 (`T_h\,\mathcal U`, `(T\times\mathcal U)`, `T_g^{-1}`, `\mathcal F/T`). Fix: write the point groups as `\mathsf T` and `\mathsf O` (as Albert et al.'s `\mathsf G` in Table D, L1077) and the transpose as `^{\top}`. Otherwise gloss both meanings in Notation.
- L40: the comment says "coordinate arrays bold", but the arrays `q`, `p` and `x` are not bold (L662-665, L680, L838, L843, L1005, L1111-1113). Also, L1067 uses plain `e_i` (a basis of a T level) next to the bold normal frame `\mathbf e_\alpha`. Fix: narrow L40 to "configuration arrays and matrices bold", and rename L1067's `e_i`, for example to `u_i`.
- Bold vectors and matrices are typed directly with no macro: `\mathbf A_s` (L688, L693-695, L699-700, L975, L1110), `\mathbf M_s` (L689, L694-695, L975, L978-979, L989, L1065, L1084, L1110), `\mathbf e_\alpha` (L661, L665, L694, L974, L1114), `\mathbf x_i` (L76, L127, L129, L1107), `\mathbf m` (L86, L93, L910, L1107), `\mathbf J` (L883) and `\mathbf n` (L983-984). Fix: add `\Amat[1]`, `\Mmat[1]` (argument form like \Pmat), `\evec`, `\xvec` and `\mvec` beside L48-50, or state in L40 that only X, P and R have macros.
- `\alpha` indexes nuclear types (L101-108, L237-244) and also normal coordinates (L661, L665, L694, L974, L1114). Fix: use another index for normal coordinates, for example `\nu` (currently unused).
- `m` denotes the masses `\mathbf m=(m_1,\ldots,m_N)` (L86, L93, L910) and the torsional label `m\in\mathbb Z` (L837-873, L1123). Fix: optional. Both rows exist in Notation (L1107, L1123), so state the two meanings there.
- `b` is the generator `(23)(45)^*` (L422) and also appears in the amino offset `b_A` (L733, L1125). Fix: rename `b_A`, for example `y_A`.
- `d` is the site separation (L630-647, L1121), `d_\Gamma` is the dimension (L521, L526, L971, L1119), and `\mathfrak d` (L1028, L1031, L1068) is the same dimension under a second symbol. Fix: write `d_{\Gamma_{\mathrm{rot}}}` at L1028, L1031 and L1068, or gloss `\mathfrak d` as Albert et al.'s symbol.
- Chemical formulas are typeset three ways: plain text "KRb + KRb" (L315), text with math subscripts "CH$_3$NH$_2$" (L415), and math roman `\mathrm{KRb}+\mathrm{KRb}` and `\mathrm K_2\mathrm{Rb}_2` (L318-320, L342). Isotopes appear as words, "carbon-12" and "nitrogen-14" (L567, L605, L917), and as a superscript, `$^{40}\mathrm K{}^{87}\mathrm{Rb}$` (L331). Fix: one form, for example `$\mathrm{KRb}$` and `$\mathrm{CH_3NH_2}$` throughout.
- L1007: `\Rmat_z(\tfrac{2\pi}{3})`, while the same angle is written `2\pi/3` at L865, L1010 and L1083. Fix: `\Rmat_z(2\pi/3)`.
- L443: `$60^\circ$` is the only angle in degrees. Torsional shifts are written in radians at L731 (`-\pi/3`) and L748 (`\pi/3`). Fix: `\pi/3`.
- L920: `\operatorname{sgn}\sigma\;\sigma\act\Psi` uses `\;`, while L259 and L712 use `\,` before `\sigma\act`. Fix: `\,`.
- L260 `\quad\forall\,\sigma` against L292 `\qquad\forall s`. Fix: `\qquad\forall\,` at both.
- L1120: Notation lists only "the species of $G_6$". It also omits other symbols used in the text: N, I_i, I_\alpha, n_\alpha, \xi, \zeta and f_\xi (Sections 1-3); G_{\mathrm{in}}, G_{\mathrm K} and G_{\mathrm{Rb}} (L324-326); v_1-v_3 (L553); D^J_{MK} and \hbar (L838); and, from the closing section, O, T, \mathrm T_{1,2}, \mathrm A, {}^1\mathrm E, {}^2\mathrm E, \mathrm T, c_2, c_3, \eta_h, \chi_{\mathrm{spin}}, \chi_{3N}, \mathfrak d, \mathfrak m_{\mathrm{nuc}}, \Gamma_{\mathrm{rot}}, \Gamma_{\mathrm{nuc}} and \mathbf n. Fix: add rows, or state the page's scope in its heading line.
- L1124: the Section 1 water symbols are the last row, while the other rows follow section order. Fix: move the row to after L1107.

## 4. LaTeX hygiene

- L451-458: Figure 5 (`[p]`) triggers "Float too large for page by 16.62926pt" (log, input line 458). Fix: `width=.96\linewidth` at L453, which removes about 24pt of height.
- L896-898: `\captionsetup{type=figure}` sits outside any group or environment, which triggers the caption warning, and the setting persists to the end of the document. Fix: `\begin{center}\captionsetup{type=figure}\includegraphics[width=\linewidth]{fig13-rigid-orientation-ball}\caption*{...}\end{center}`. Keep `\caption*`, because build.py counts `\caption{`.
- L1090: `\FloatBarrier` is redundant twice over. `\clearpage` at L1089 already flushes all floats, and placeins[section] adds a barrier at `\section*{Notation}` (L1094). Fix: delete L1090.
- L61, L888, L1089, L1129: manual `\clearpage` after the guide plate and before the closing section, the Notation and the bibliography. Fix: none needed if each unnumbered part should open a page. State that rule in the L26 comment block.
- L903: the whole closing section is wrapped in `{\small ...}`, with display skips cut to 4pt/2pt and `\parskip` 2pt. This is an inline page-fitting override. Fix: define it once in the preamble as a named environment, or drop it and let the text flow (the brief says layout is not final).
- L904: \closingtable carries `\\[3pt]` (outside any tabular) and `\par\medskip`. Fix: see section 2, C1.
- L905, L909, L917, L923, L935, L941, L954, L960, L973, L992, L1002, L1054: there are 12 hand-built run-in paragraphs (`\noindent\textbf{...}\quad`). Fix: `\titleformat{\paragraph}[runin]{\normalfont\bfseries}{}{0pt}{}` with `\titlespacing*{\paragraph}{0pt}{.5ex}{1em}`, then `\paragraph{...}`. For L905, set `\parindent` to 0pt in the closing environment.
- L482/483, L726/728, L862/863: two itemize lists run back to back, which doubles \topsep. Fix: merge each pair into one list (or note that the split is intended).
- L948-949: `\mathrm` ends L948 and its argument `E` starts L949. Fix: keep `\mathrm E` on one line.
- L890 against L892: the banner says "(unnumbered, three pages)" and then "flow over four pages". The build puts the closing section on pp. 20-23, which is four pages. Fix: "(unnumbered, four pages)", or drop the page counts.
- L892: this banner line runs far past the roughly 100-column width of every other banner line. Fix: rewrap.
- Comments: there are no references to deleted files (docs/caption-notes.md or any other .md) and no history narration. compute/make_data.py and checks/verify_methane.py (L892-893) exist. The L10 comment is covered in section 1. No action.
- Unused labels: sec:feasible-regimes (L349), sec:localized-states (L610) and sec:position-coordinates (L659) are never referenced. fig:1-7, fig:9, fig:10 and fig:11 are never `\ref`ed, but build.py check_numbering reads them, so keep them. Fix: use sec:localized-states at L970 (see the hard-coded references item below), and either reference or delete the other two sec labels.
- L349 and L659: two label names are stale. sec:feasible-regimes labels "Symmetry Groups and Versions", and sec:position-coordinates labels "Position Representation". Fix: sec:groups-versions and sec:position (neither is referenced, so renaming breaks nothing).
- Hard-coded references: L902 "Figure~11" should be `Figure~\ref{fig:11}`. L970 "Section~9" should be `Section~\ref{sec:localized-states}`. L1009 "Section~12" should be `Section~\ref{sec:momentum}`. L1011 "Figure~12" should be `Figure~\ref{fig:12}`. L1124 "(figure 1)" should be `(Figure~\ref{fig:1})`.
- L905 and the run-in numbers at L909 ("1."), L917 ("2--3."), L923 ("4--5."), L935, L941, L954, L960, L973, L992 and L1002 hard-code section numbers. Fix: `\ref`. This needs new labels on subsections 1, 2, 3, 4 and 11 (L69, L163, L229, L311, L779).
- L925, L951, L954, L958, L1013, L1060, L1062, L1064 and L1082: Tables A-D are cited by hard-coded letters because \closingtable has no counter. Fix: see section 2, C1 (counter, `\label`, `Table~\ref{...}`).
- Missing tie before `\cite` at L105, L106 and L408 (line break before `\cite`), at L333 and L808 (`et al.\ \cite`), and at L441-442 (`et al.\` then a newline). Ties are used at L442 (Kr\k{e}glewski), L444, L739, L796, L810 and L1055. Fix: `~\cite` everywhere.
- L796: "Section 1.3" in the Hatcher footnote has no tie, while every other "Section" is tied. Fix: `Section~1.3`.
- Status, no action: all seven `\ref` uses are tied (L446 twice, L584, L639, L755, L872, L875). Quote marks: the only pair is TeX-style ``...'' at L1069.
- L684: the display ends with "," but the next `\item` (L687) starts a new sentence. Fix: ".".
- L573: the display ends with ",", and the sentence continues as a new `\item` that starts with lowercase "so" (L576). Fix: end the display with "." (the author decides the wording of L576), or keep the continuation inside one item.
- L481-482 and L807-808: "the example that closes the paper" is a cross-reference in words with no link. Fix: optional. Add `\label{sec:example}` after L895 and wrap the words in `\hyperref[sec:example]{...}`, which leaves the wording unchanged.

## 5. Heading format

Status (complies): L22 renders parts as "I. Configurations and States" through "IV. ...". L24 renders "1. Configuration Space" through "12. Momentum Representation". L23 renders the unnumbered titles bare: "Example and Recovery of the Rigid Case", "Notation" and "References". No heading has a trailing period. The run-in heads at L909-1002 and L1054 have no trailing period. The period in the table titles "Table A." etc. is a caption-label separator; section 2 covers it (labelsep=period).

- L909, L917, L923, L935, L941, L954, L960, L973, L992, L1002 and L1054: the run-in heads are in sentence case ("1. Configuration space"), while every numbered and unnumbered heading is in Title Case ("1. Configuration Space", L69). Fix: Title Case in the run-in heads, or record sentence case as the run-in convention in the preamble comment at L21. The heads are hand-built; see the run-in paragraph item in section 4.
- L36: the title "...Molecules and their Complexes" lowercases "their", while the headings use Title Case. The title also contains a manual `\\` break. Fix: "Their", if the title follows the headings' Title Case.
