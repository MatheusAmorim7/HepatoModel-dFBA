# Zonação enzimática no PhysiCell nativo (etapa 1): resultados

Exploratório; sem tumor; NÃO é calibração nem validação. Escopo registrado (`postabsorptive-native-zonation-2026-11-21`), protocolo e previsões em `PROTOCOL_Z.md` (escritos antes de rodar). Configuração: pós-absortivo, Q=0,86 L/min, dt=0,1 min, estoque local de sangue, 30 min, 3000 hospedeiros, transporte serial periportal→central (Fase C).

## O que foi implementado
Cópia do nativo (`native_z`, patch de 1 linha: peso de tarefa por hospedeiro). Para cada hospedeiro, pelo segmento do seu voxel (0=PP…7=central): tarefa de exportação de glicose 2,14·m (e teto de glicogenólise·m), tetos de PEPCK/FBP1/G6Pase = κ·m·J_ref (κ=1,5), com m = abundância linear PP-alta (s=3 ou 10), média 1 sobre os hospedeiros. Total do órgão preservado. Zonação do fluxo gliconeogênico é **prescrita**; emergem os gradientes de substrato, extração e viabilidade.

## Placar das previsões
| | Resultado |
|---|---|
| P1 Z0 = cL3 por agente ≤1e-6 | **FALHOU como escrita**: 0,27 por hospedeiro (a partir do passo 19; 225 hospedeiros). Soma por voxel idêntica a 3e-14 em todos os passos: só muda a divisão degenerada de estoque limitante entre hospedeiros do mesmo voxel (os tetos diferem entre os dois, o que altera o vértice do LP). Balanços ≤8,4e-13 pmol. |
| P2 Z3c (só tetos) rejeitado no passo 0 em segmento ≥4 | **CONFIRMADA**: LP local inviável, voxel 972 (segmento 6). Sem tarefa por posição, tetos zonais são inviáveis: a divisão do trabalho tem de ser tarefa. |
| P3 Z3, Z10 viáveis 30 min, totais do órgão ≤3% de Z0 | **CONFIRMADA**: viáveis; totais de lactato, glicerol e O₂ idênticos a Z0 (0,0%); balanços ≤8,3e-13. |
| P4 PEPCK monotônico, razão seg0/seg7 dentro de 10% de 2,56 | **PARCIAL**: monotônico (Z3, Z10); razão 2,86 (Z3; +12%) e 9,63 (Z10). O excesso vem do glicerol (G1, teto fixo e extraído no PP), que alivia o PEPCK no PP. |
| P5 queda de lactato PC/PP mais forte em Z3 que em Z0 | **FALHOU**: −17,5% (Z3), −15,4% (Z10) contra −21,4% (Z0). A extração total é igual, mas fica concentrada no PP e o segmento 0 já está mais esgotado; a razão vs segmento 0 achata. Métrica mal escolhida por mim. |
| P6 glicerol: extração ≥0,9 mantida; fração PP maior em Z3 | **PARCIAL**: extração 0,946 em todos; fração nos segmentos 0–2 = 0,672 em todos (teto de glicerol não zonado). |
| P7 estado estacionário ≤6 min, deriva <1% | **CONFIRMADA**: estacionário em ≈2–3 min; deriva 0,0%. |

## Resultados centrais (Z3 vs Z0)
- Captação de lactato por hospedeiro: Z0 1,35→1,40 (plano); Z3 1,75 (PP)→0,52 (PC); O₂: Z3 2,58→1,75. Glicose (G6Pase) 2,50→0,98 por hospedeiro.
- Totais do órgão iguais: a captação é fixada pela tarefa e pelo suprimento (fluxo-limitado); a zonação redistribui onde, não quanto.
- Glicerol esgota no PP (segmento 0→6: queda de 99%); amônia, alanina e palmitato quase planos (−6% a −3%).
- Tudo estável por 30 min em Q=0,86; Q=0,43 continua sem sustentar (Fase C).

## Limites
Perfil linear e s, κ são hipóteses (só o 3× de PK em rato tem fonte); ureia/GLUL/GK/PK não zonados; hospedeiro-a-hospedeiro dentro do voxel é degenerado (usar agregados por voxel/segmento); sem validação independente; revisão independente pendente (dispatcher indisponível). Figuras em `figures/`.

---
# ADENDO E: fluxo gliconeogênico EMERGENTE (2026-11-22; aprovação `postabsorptive-native-zonation-emergent-2026-11-22`)
Sem tarefa de glicose: cada hospedeiro maximiza a exportação de glicose (recompensa λ=200, binário `native_ze`) sujeito a capacidade (tetos de PEPCK/FBP1/G6Pase = κ·m·J_ref; glicogenólise ·m), suprimento, ureia 0,19 e ATPM 5. Entrada = abundância enzimática (hipótese linear s=3); saída (fluxo) = previsão. Q=0,86 L/min, dt=0,1, 30 min; todos os 4 braços viáveis, balanços ≤8,7e-13 pmol, sem rejeição.

| Braço | Glicose média/hospedeiro | PEPCK seg0/seg7 (capacidade) | Extração de lactato | vs Z0: lactato / O₂ |
|---|---|---|---|---|
| E0 uniforme, κ=1 | 2,119 | 1,00 (1,00) | 0,61 | −3,0% / −1,2% |
| E3 s=3, κ=1 | 2,116 | 2,78 (2,56) | 0,60 | −3,5% / −1,4% |
| E0 uniforme, κ=1,5 | 2,458 | 1,03 (1,00) | 0,91 | +46% / +18% |
| E3 s=3, κ=1,5 | 2,453 | 2,54 (2,56) | 0,91 | +45% / +18% |

## Placar (Q1–Q6)
- **Q1 (E0_k10 = referência a ±0,02 e totais ≤2%)**: **FALHOU por pouco**: 2,119 (−1,0%), lactato −3,0%; os hospedeiros chegam a 95–99% do teto de FBP1 (não 100%), provável limite local de troca tecido–sangue/estoque por passo (não verificado).
- **Q2 (E3_k10 viável, 2,14±3%, tudo no teto)**: **PARCIAL**: viável e exportação 2,116 (−1,1%); mas PEPCK fica abaixo do teto no pericentral (0,92 no seg 7), então o gradiente de fluxo (2,78) é um pouco mais íngreme que o de capacidade (2,56): efeito de suprimento.
- **Q3 (E3_k15: +10–40% de exportação, extração >0,8, PC abaixo do teto)**: **PARCIAL**: +14,6% e extração 0,91 confirmados; mas o déficit de capacidade está no periportal (FBP1 0,90–0,91 do teto), não no pericentral; razão do fluxo 2,54 ≈ capacidade (não mais íngreme).
- **Q4 (E0_k15: PP>PC por suprimento, razão ≥1,2)**: **FALHOU**: razão 1,03; o déficit aparece no meio do lóbulo (segmentos 4–6, PEPCK 0,91–0,93 do teto, FBP1 0,84–0,86) com lactato no sangue a 5–6% da entrada, e recupera no centro.
- **Q5**: **CONFIRMADA** (deriva 0,0%; estado estacionário em 4 min nos braços κ=1,5, imediato nos κ=1). **Q6**: **CONFIRMADA**.

## Leitura científica
1. Com a recompensa, o gradiente de fluxo **segue a capacidade enzimática** (E3_k10 ≈ prescrito Z3 nas curvas): o suprimento só o modifica em ≤10%. Logo "emergir" move a hipótese de "fluxo" para "abundância enzimática" (mensurável por proteômica), mas não gera zonação sem ela: com enzimas uniformes o fluxo fica plano (razão 1,0–1,03).
2. O que emerge de verdade é onde o suprimento limita: com capacidade excedente (κ=1,5) a extração de lactato sobe a 0,91 (o limite de oferta a Q=0,86) e a produção de glicose sobe só 14,6% (vs +50% de capacidade); com Q menor isso cortaria a produção (Fase C: a 0,43 não sustenta).
3. O glicerol (extração 0,946, esgotado antes do segmento 5) e os totais de glicerol não mudam entre braços.
Limites: perfil linear, s e κ são hipóteses; λ=200 corresponde a um regime de "máxima produção" (não calibrado); glicogenólise no teto em todos os segmentos (a cap ·m também é hipótese); ureia/GLUL/GK/PK não zonados; revisão independente pendente. Figuras fig6–fig8. Scores em `SCORES_ZE.json`.
