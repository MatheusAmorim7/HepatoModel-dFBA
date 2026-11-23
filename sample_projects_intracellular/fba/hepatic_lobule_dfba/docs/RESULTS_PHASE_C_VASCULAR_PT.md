# Fase C: transporte multi-campo em sinusoide serial, estoques finitos e passo metabolico (relatorio)

Escopo aprovado: id postabsorptive-native-phase-c-2026-11-20 (autorizacao: "aprovo, vamos seguir"; registro oficial via state.py; gate PASSOU). Exploratorio; nao e calibracao nem validacao; sem tumor, zonacao enzimatica, ciclo/morte/movimento.
Protocolo e previsoes: PROTOCOL_C.md (com 5 adendos pos-hoc declarados). Sem sessao formal. Todos os casos planejados foram concluidos (cE2 incluido). Revisao independente (phase_c_review) pendente (dispatcher indisponivel nesta sessao).

## O que foi construido
- src/vascular.h: operador de transporte de UM campo (Euler implicito; tecido 50x50 em banda LAPACK + complemento de Schur sobre 48 nos sanguineos em serie, 6 setores x 8 segmentos do mapa candidato; segmento 0 = anel externo periportal, 7 = centro). Verificado contra numpy (~1e-14, 3 regimes) e contra solucao densa (1,8e-15).
- native_c (derivado de native_b0b por make_native_c.py): transporte para os 25 campos com estoques finitos, entrada/saida, balanco certificado por campo, vazao/volume/passo por caso, tampao linear de O2 (hemoglobina; beta=107,45, hipotese exploratoria).
- native_c2 (derivado por make_native_c2.py): estoque local de sangue (LBP): o LP ve o estoque do tecido + a parcela do voxel no conteudo do seu no; captacao alem do tecido e debitada do no (conserva massa).
- Parametros VARRIDOS/HERDADOS, nao adotados: vazao 0,43/0,86/10 L/min (18,75 nL/min por L/min), volume sanguineo 2 nL (1 nL em cD), ge=gw=3,504 nL/min herdados, entradas de lactato/palmitato/glicerol/aminoacidos/amonia citadas (PROTOCOL_C.md).

## Resultados por criterio
| Criterio | Resultado |
|---|---|
| 1 Balanco por campo <= 1e-9 pmol, sem negativos/reparo | PASSOU em todos os casos: erro maximo 1,5e-12 pmol, residuo de transporte maximo 3,2e-13; transacoes rejeitadas sao atomicas (sem reparo) |
| 2 Refinamento de passo muda < 1% ao dividir por 2 | FALHOU em parte: convergencia de PRIMEIRA ORDEM (erro cai pela metade a cada divisao): sem LBP a 2 min, 0,04->0,02: O2 1,6%, lactato 1,0%, glicerol 6,9%; 0,02->0,01: O2 0,8%, lactato 0,5%, amonia 0,17%, aminoacidos < 0,2%, glicerol 3,4%. Com LBP a 5 min, 0,1->0,05: O2 4,2%, lactato 5,8%, glicerol 12,8%, amonia 1,3%, palmitato 1,5%, alanina 0,6%. O estado estacionario do glicerol difere 3,7% entre dt=0,1 e 0,05 (1,129 vs 1,172 pmol/min); dt=0,05 coincide com dt=0,04 sem LBP (1,171) |
| 3 Estado estacionario = estatico com oferta suficiente | PASSOU: cA2 (dt=0,02) e cL1 (dt=0,1, LBP): desvio 0,0 do B0b em todas as 225 reacoes e passos; cL5 (LBP) e cB_q086 identicos (0,0) e o LBP nunca foi acionado (retirada do no = 0) |
| 4 Escassez classificada, sem reparo | PASSOU: O2 so dissolvido rejeitado em 0,08 min; Q=0,43 com LBP (cL4) rejeitado aos 3,1 min (escassez real: glicerol 1,38->0,61 pmol/min contra entrada 0,60, lactato e O2 compensando); rejeicoes em dt maior (cE_q043_long dt=0,04 aos 0,24 min; cL2 dt=1,0; cL6 dt=0,2) sao ARTEFATO de passo, nao escassez |
| 5 Perfil portal->central e sensibilidade | PASSOU em parte: a Q=0,86 (cL3, 30 min) o lactato cai 21% do segmento 0 ao 7, glicerol 90%, amonia 6%, palmitato 6%, O2 4%, alanina 2,7%; sensibilidade a vazao 0,43 vs 0,86 mostrada (0,43 rejeita, 0,86 se sustenta); sensibilidade ao volume sanguineo (2 vs 1 nL) INCONCLUSIVA: ambos ainda no estado estatico aos 2 min |

## Placar das previsoes (revisadas e posteriores)
| Previsao | Resultado |
|---|---|
| P1' cA2 reproduz B0b | CONFIRMADA (0,0) |
| P2' O2 dissolvido rejeitado em ate 0,5 min | CONFIRMADA (0,08 min; extracao de O2 impossivel de 10,9) |
| P3' Q=0,43 rejeitado em 0,3-2 min; Q=0,86 viavel com glicerol abaixo do teto e lactato 0,5-0,95 | PARCIAL/FALHOU: Q=0,43 completou 2 min (erro de escala de tempo); Q=0,86 viavel e lactato 0,61 (ok); glicerol ainda no teto em 2 min (cai a 1,17 entre 2 e 4 min) |
| P4' aminoacidos/NH3 < 1% (0,02->0,01); O2/lactato/glicerol >= 1% (0,04->0,02) | CONFIRMADA |
| P5' lactato cai >= 10% PP->PC; V=1 nL aproxima mais rapido | PARCIAL: lactato -13% (2 min) e -21% (30 min) confirmado; volume inconclusivo |
| P6 Q=0,43 sem LBP rejeita entre 6 e 20 min | FALHOU: cE2 (Q=0,43, dt=0,02, sem LBP, relancado apos o reinicio do contêiner) foi REJEITADO aos 3,56 min (178 passos aceitos), sem artefato de passo; o glicerol fica no limite da entrada (0,596 contra 0,597 pmol/min), o lactato a 133% da entrada (puxando estoque do tecido) e o LP local fica inviavel quando o estoque acaba. Concorda com cL4 (LBP, dt=0,1: 3,1 min) |
| P7 Q=0,86: glicerol abaixo do teto em ate 5 min, sem rejeicao ate 10 min | CONFIRMADA (1,38->1,17 pmol/min entre 2 e 4 min, extracao 0,98 da entrada; estado estacionario em ~4 min) |
| P8 afastamento do estatico (> 0,1) antes de qualquer rejeicao | CONFIRMADA (0,40 nos casos com LBP; 1,02 sem LBP a dt=0,04, onde o artefato de passo soma) |
| PL1 cL1 e cL2 reproduzem o estatico | PARCIAL: cL1 (dt=0,1) confirmada; cL2 (dt=1,0) REJEITADO no primeiro passo |
| PL2 cL5 = cB_q086 | CONFIRMADA (fluxos e perfis identicos) |
| PL3 Q=0,86 30 min: glicerol abaixo do teto em ate 10 min e <= 1,19; lactato > 8,02; sem rejeicao | CONFIRMADA (glicerol 1,129, lactato 8,28, estacionario de 6 a 30 min) |
| PL4 Q=0,43 com LBP rejeita entre 5 e 30 min | FALHOU: rejeitado aos 3,1 min (mais cedo; o glicerol acaba rapido) |
| PL5 refinamento com LBP < 1% por divisao a 5 min | FALHOU (ver criterio 2) |
| PL6 reajuste de arredondamento < 1e-12 pmol | CONFIRMADA (0 em todos os casos) |

## Achados
1. A hemoglobina e indispensavel: so com O2 dissolvido a entrada cobre 5-9% do consumo; com Hb a extracao de O2 e 10-20%.
2. O passo e limitado por separacao de operadores: o LP so enxerga o estoque do voxel no inicio do passo. Sem LBP o passo util e <~0,04 min; com LBP, ~0,1 min; a dt=0,2 e 1,0 o LBP nao basta. Chegar a minutos exige acoplamento monolitico (consumo como termo-fonte na solucao implicita do transporte).
3. A precisao quantitativa pede dt <= 0,05 (convergencia de primeira ordem); dt=0,1 serve a analises qualitativas e a horizontes longos.
4. O cenario pos-absortivo (S7, glicose 2,14, ureia 0,19) se sustenta a 0,86 L/min, com glicerol limitado pela oferta (extracao 0,95-0,98) e lactato com extracao ~0,63, e NAO se sustenta a 0,43 L/min: a escassez aparece aos ~3-3,5 min de forma robusta (3,1 min com LBP e dt=0,1; 3,56 min sem LBP e dt=0,02), bem antes do previsto (6-20 min).
5. O estado estacionario e atingido em ~4-6 min e permanece estavel por 30 min (cL3).

## Limites
Condutancias herdadas (acoplamento sangue-tecido quase instantaneo: regime limitado por fluxo); volume sanguineo, beta de O2 e vazao sao hipoteses varridas; glicerol/lactato/palmitato portal = arterial assumido; o horizonte alcancado foi 30 min (nao 60) e o passo util <= 0,1 min (nao 1 min); 3000 hospedeiros com 1/3 de variabilidade de limites zonais so nas cinco enzimas herdadas; sem tumor, sem ciclo/morte/movimento; sem revisao independente.
Reproducao: make -f .../Makefile (native_c, native_c2); make_inputs_c*.py; run_cases*.sh; analyze_case.py <caso>; compare_refinement*.py.
