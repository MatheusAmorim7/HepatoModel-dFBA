# Calibração do ATP de manutenção contra o consumo hepático de O₂ (2026-11-22)

**Calibração declarada de UM parâmetro (CULTURE_ATPM 5 → 21,6) contra UM alvo** (O₂ hepático 42,9 mL/min = Elia/Wang 200 kcal/kg/dia × 1,5 kg ÷ 4,86 kcal/L; valores de resumos de busca, originais não lidos). Tudo o mais é validação não usada no ajuste. Exploratório. Protocolo e previsões: `PROTOCOL_ATP.md` (adendos A e B escritos antes das rodadas correspondentes). Aprovações: `postabsorptive-atp-calibration-2026-11-22` e `...-amend-...`.

## Cronologia honesta
1. Scan estático (scan_atpm.py, L1, cenário de referência): ATPM=21,6 dá 42,94 mL/min; o teto de O₂ antigo (−6) truncava em 43,1 mL/min e foi levantado para −20 (mudança declarada). ATPM máximo viável ≈23 com o teto antigo.
2. 1ª tentativa no nativo: falhou (3 braços, rc=1) por uma verificação fixa `ATPM==5` no código; criei cópia com a verificação parametrizada (`native_zatp`, adendo A).
3. 2ª tentativa (dt=0,1, 30 min): 3 braços rejeitados no passo 0 (consumo de O₂ de um passo > O₂ disponível; limite de passo da Fase C, 2,6× pior). Preservado em `runs/failed_dt01_step0`. Adendo B: dt=0,02, horizonte 6 min.
4. 3ª tentativa (dt=0,02, 6 min): os 3 braços viáveis, balanços ≤6,7e-14 pmol, sem rejeição.

## Placar (R1–R7)
| | Resultado |
|---|---|
| R1 O₂ do órgão 42,9 ±5% | **CONFIRMADA**: 43,04 mL/min (estático 42,94). Ajuste, não validação. |
| R2 A1 (só 0,124 mM de palmitato) rejeitado antes de 10 min | **FALHOU**: A1 viável 6 min; a extração de palmitato é 0,76 e basta. O palmitato "carbono-equivalente" (A2) não mudou nada (A1 = A2 em todos os totais). |
| R3 queda do conteúdo de O₂ PP→PC ≥20% (era ~4%) | **CONFIRMADA por pouco**: 20,5% (A2), 20,1% (A3). |
| R4 palmitato 81±20 µmol/min, abaixo da faixa de Iozzo | **CONFIRMADA**: 81,5 µmol/min (abaixo de 120–240; unidades não reverificadas). |
| R5 lactato, glicerol, glicose ≤5% vs ATPM=5 | **PARCIAL**: lactato −0,7%, glicose 2,14; glicerol +5,3% (marginalmente fora). |
| R6 A3: O₂ por hospedeiro PP/PC ≥1,4, sem rejeição | **FALHOU na razão**: 1,14 (6,25 → 5,47); sem rejeição. A carga de ATP comum domina; a tarefa zonada de glicose contribui pouco. |
| R7 balanços ≤1e-9 | **CONFIRMADA** (≤6,7e-14). |

## Validação não ajustada (não usada no ajuste)
- **Gradiente de O₂ no sinusoide:** conteúdo cai de 0,93 para 0,74 da entrada; PO₂ estimada* 56 → 37 mmHg, contra 60–65 → 30–35 mmHg na literatura ([Thieme](https://www.thieme-connect.de/products/ejournals/html/10.1055/s-0041-1742279); [NCBI NBK53069](https://www.ncbi.nlm.nih.gov/books/NBK53069/); parte em rato). *Conversão por curva de Hill externa ao modelo (P50 26,8, n 2,7, saturação de entrada 0,95): aproximação. Antes da calibração o gradiente era ~57 → 51 mmHg (quase plano).
- **Extração hepática de O₂:** 0,26 vs ~0,37–0,40 relatado (resumos: Deranged Physiology, Turku PET; valor secundário). Modelo abaixo, coerente com um VO₂ alvo menor que o de Lutz (90 mL/min).
- **Palmitato:** 81 µmol/min, provavelmente abaixo do dado de Iozzo (não conclusivo).
- **Lactato:** extração 0,62 (não calibrada).

## Leitura científica
1. Com um único parâmetro (ATPM) calibrado ao consumo de O₂, o gradiente de O₂ no sinusoide passa de ~4% para ~20% e cai na faixa da literatura sem ter sido ajustado a ele: é o resultado de validação mais útil até aqui, embora aproximado (Hill externa).
2. Não surgiu hipóxia: O₂ por hospedeiro permanece ~6 em todos os segmentos (sem limitação por O₂) e a zonação de glicose quase não muda o consumo de O₂ (razão 1,14). O consumo de O₂ é dominado pelo ATP comum.
3. Limite estrutural: o alvo alto de literatura (≈90 mL/min) exigiria ATPM ≈50, além do que os tetos de combustíveis atuais sustentam no estático (acima de ATPM≈40 o L1 passa a oxidar mais lactato; palmitato encosta em 0,4). Não testado no nativo.
## Limites
Alvo único com incerteza de fator ~2 entre fontes; ATPM é um parâmetro agregado (síntese de proteína, bombas iônicas, VLDL não estão separados); alvo e validações vêm de resumos de busca (conferir originais); horizonte 6 min; dt=0,02; zonação prescrita (s=3) herdada; revisão independente pendente.
Arquivos: SCORES_ATP.json, SCAN_ATPM.json, figures/fig_atp1_gradiente_O2.png, fig_atp2_vo2.png.
