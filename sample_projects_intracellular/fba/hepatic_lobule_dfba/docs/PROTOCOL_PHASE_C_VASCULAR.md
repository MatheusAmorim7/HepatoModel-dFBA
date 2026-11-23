# Protocolo Fase C: transporte serial multi-campo (registrado antes de gerar entradas e rodar)

Escopo aprovado: id postabsorptive-native-phase-c-2026-11-20 (autorizacao: "aprovo, vamos seguir"). Exploratorio; nao e calibracao. Sem tumor, zonacao enzimatica, ciclo/morte/movimento.

## Implementacao (ja construida e verificada)
- src/vascular.h: operador de transporte de UM campo (Euler implicito; tecido 50x50 em banda LAPACK + complemento de Schur sobre 48 nos sanguineos em serie, 6 setores x 8 segmentos do mapa candidato de
  lobule_vascular_geometry_bridge_2026_11_17; segmento 0 = anel externo periportal, 7 = centro). Verificado contra numpy (diferenca ~1e-14 em 3 regimes) e contra solucao densa (1,8e-15); certificado de equacao <= 1e-10 e
  de massa <= 1e-9; sem clipping (negativo reprova).
- src/cadence/native_c.cpp (derivado de native_b0b.cpp por make_native_c.py, com substituicoes verificadas): transporte para os 25 campos com estoques finitos, vazao e volume por caso, balanco por campo, intervalos de fenotipo/mecanica por caso.
## Parametros (varridos ou herdados; NAO adotados como fisiologicos)
- vazao plasmatica do dominio = Q(L/min) x 18,75 nL/min por L/min (3000 hospedeiros representando 1,6e11 hepatocitos: 0,43 L/min -> 8,06 nL/min); casos Q = 10 (controle nao limitante), 0,86, 0,43;
- volume sanguineo sinusoidal total 2 nL (10% do volume geometrico de 20 nL; tempo de transito V/Q = 0,25 min a 0,43 L/min... 0,12 min a 0,86), sensibilidade 1 nL;
- condutancias herdadas do banho: ge = gw = 3,504 nL/min por voxel (nao medidas; acoplamento tecido-sangue quase instantaneo, regime limitado por fluxo);
- concentracoes de entrada (mM), misturando portal e arteria: aminoacidos de Neis (Tabela 2), lactato 0,821 (Owen, arterial = portal assumido), palmitato 0,124 (CHS, idem), glicerol 0,074 (coorte nao obesa, valor de jejum),
  amonia 0,239 (derivada: entrada ~103 umol/min de Levitt / 0,43 L/min; inferencia), O2 herdado do controle; campos de exportacao e sem fonte = 0; tecido inicial = concentracao de entrada.
- passo metabolico dt em minutos; estoques finitos reais (sem 1,0 pmol); cache de LP morno.
## Casos e criterios (fixos agora)
A (controle de suprimento nao limitante): Q=10, dt=0,1, 20 passos. Criterio 3: fluxos aplicados no ultimo passo = B0b (cenario estatico) <= 1e-6 por agente em todas as 225 reacoes; sem rejeicoes.
B (suprimento fisiologico): Q=0,86 e Q=0,43, dt=0,1, ate 60 passos (6 min): reportar viabilidade, rejeicoes atomicas (classificadas), uso de glicerol/lactato/alanina vs oferta, extracao por campo, perfil portal->central por campo.
C (refinamento): Q=0,86, dt in {0,2; 0,1; 0,05} no mesmo horizonte de 2 min; mudanca relativa de captacoes totais e estoques ao dividir o passo por 2 < 1% por campo, ou reportar a falha.
D (sensibilidade): volume sanguineo 1 vs 2 nL a Q=0,86.
Criterios transversais: balanco por campo <= 1e-9 pmol por passo; sem negativos; sem reparo; passos aceitos so por transacao atomica; arquivos selados e de entrada inalterados por hash.
## Previsoes (escritas antes de rodar)
P1 (alta): caso A reproduz o B0b a <= 1e-6 e fecha todos os balancos.
P2 (media): em Q=0,43 o LP local fica inviavel (rejeicao atomica classificada) em tempo simulado entre 0,3 e 3 min (estoque inicial de lactato e de glicerol se esgota; a oferta nao cobre o requisito de lactato, extracao estatica 1,21).
P3 (media): em Q=0,86 o cenario permanece viavel por 6 min, com glicerol abaixo do teto (limitado pela oferta, extracao <= 1) e lactato com extracao entre 0,5 e 0,95.
P4 (media): refinamento: aminoacidos e NH3 mudam < 1% entre dt 0,1 e 0,05; lactato e glicerol mudam >= 1% entre dt 0,2 e 0,1 (divisao do passo reduz o artefato de separacao de operadores).
P5 (media): em Q=0,86 o lactato cai >= 10% do segmento 0 ao 7 em estado quase estacionario; o volume sanguineo 1 nL aproxima o estado estacionario mais rapido que 2 nL.
Se alguma previsao falhar, reportar a falha; sem ajustar parametros para obter viabilidade.

## ADENDO pos-hoc (escrito DEPOIS de ver a execucao rejeitada de cA_nonlimiting a dt=0,1; esse resultado e preservado em runs/cA_nonlimiting)
Fato observado: o caso A (suprimento nao limitante, dt=0,1) foi rejeitado atomicamente no primeiro passo (voxel 69, 3 hospedeiros; LP local inviavel; rejection.json). Causa verificada: o O2 dissolvido no voxel (0,0787 mM x 0,008 nL
= 6,3e-4 pmol) e menor que o consumo de UM passo de 3 hospedeiros (1,4e-3 pmol a dt=0,1 min); o LP so enxerga o estoque do voxel no inicio do passo, nao a reposicao do sangue durante o passo (limite de separacao de operadores).
Alem disso, o O2 dissolvido sozinho entregaria 5-9% do consumo do dominio (0,63-1,27 pmol/min contra 13,85 pmol/min); a hemoglobina e o que transporta o O2 (conteudo ~8,9 mM x 0,95 de saturacao, capacidade padrao assumida: 15 g/dL x 1,34 mL/g).
Mudancas declaradas (pos-hoc, nao confirmatorias):
- tampao linear de O2 no sangue: armazenamento e fluxo dos nos de O2 escalados por (1+beta) com beta=107,45 = (c_diss + 8,9 x 0,95)/c_diss - 1 (linearizacao na entrada; NAO e cinetica de Hb); controle beta=0 (so dissolvido).
  O tecido continua com O2 dissolvido. Parametro exploratorio, varrido (0 e 107,45), nao adotado como fisiologico.
- passos dt em {0,04; 0,02; 0,01} min (a razao consumo/estoque do voxel mais denso e 0,88; 0,44; 0,22), horizonte 0,5-2 min simulados; a regra consumo/estoque limita dt, nao a vontade de simular 60 min.
Casos revisados: cA2 (Q=10, beta=107, dt=0,02, 25 passos; controle de suprimento); cA3 (Q=0,86, beta=0, dt=0,02, 25 passos; O2 dissolvido); cB_q086 e cB_q043 (beta=107, dt=0,02, 100 passos = 2 min); cC refinamento (Q=0,86; dt=0,04 x 50 passos e
dt=0,01 x 200 passos, 2 min); cD (V=1 nL, Q=0,86, dt=0,02, 100 passos). Mesmos criterios transversais.
Previsoes revisadas (escritas antes de rodar os casos novos):
P1' (alta): cA2 reproduz o B0b a <= 1e-6 em todos os passos e fecha os balancos.
P2' (media-alta): cA3 e rejeitado atomicamente por escassez de O2 em ate 0,5 min simulados.
P3' (media): em Q=0,43 (com Hb) a rede e rejeitada por escassez (lactato/glicerol) entre 0,3 e 2 min; em Q=0,86 permanece viavel por 2 min com glicerol abaixo do teto e extracao de lactato entre 0,5 e 0,95.
P4' (media): aminoacidos e NH3 mudam < 1% entre dt=0,02 e 0,01; O2/lactato/glicerol mudam >= 1% entre dt=0,04 e 0,02.
P5' (media): em Q=0,86 o lactato cai >= 10% do segmento 0 ao 7; V=1 nL aproxima o estado quase estacionario mais rapido que 2 nL.

## ADENDO 2 (escrito DEPOIS dos resultados dos casos de 2 min; placar honesto das previsoes revisadas)
- P1' CONFIRMADA: cA2 = B0b com diferenca 0,0 (primeiro e ultimo passo), balanco por campo 2,2e-13 pmol.
- P2' CONFIRMADA: cA3 (O2 so dissolvido) rejeitado atomicamente com 4 passos aceitos (0,08 min); extracao de O2 impossivel de 10,9.
- P3' FALHOU para Q=0,43: completou 100 passos (2 min) sem rejeicao (previ rejeicao entre 0,3 e 2 min). Erro de escala de tempo: o deficit de lactato (~1,4 pmol/min) e pequeno frente ao estoque de tecido (~16 pmol).
  Para Q=0,86: viavel por 2 min (confirmada), extracao de lactato 0,61 (dentro de 0,5-0,95, confirmada); glicerol NAO ficou abaixo do teto em 2 min (consumo 1,38 > entrada 1,19 pmol/min, coberto por estoque): falhou nesse ponto.
- cC_q086_dt004 NAO rodou por erro do gerador (dt_mechanics=0,1 nao e multiplo inteiro de 0,04; Cadence_integers, rc=1); corrigido com dt_mechanics=0,12. Preservado em runs/cC_q086_dt004_stderr.txt.
Casos longos (pos-hoc, definidos depois de ver os resultados): dt=0,04, Q=0,86 por 10 min (250 passos) e Q=0,43 por 20 min (500 passos); refazer cC_q086_dt004 (dt=0,04, 50 passos).
Previsoes (escritas antes de rodar):
P6 (media): Q=0,43 sofre rejeicao atomica por escassez (lactato ou glicerol) entre 6 e 20 min simulados (estimativa: estoque de ~18 pmol de lactato / deficit ~1,4 pmol/min ~ 13 min; o glicerol pode esgotar antes).
P7 (media): em Q=0,86 o glicerol cai abaixo do teto do cenario (captacao < 1,38 pmol/min) em ate 5 min, sem rejeicao ate 10 min (o LP substitui por lactato dentro da oferta).
P8 (media): ao esgotar o estoque, o fluxo aplicado se afasta do cenario estatico (diferenca maxima > 0,1 unidades) antes de qualquer rejeicao.

## ADENDO 3 (escrito DEPOIS da rejeicao de cE_q043_long a dt=0,04)
cE_q043_long (Q=0,43, dt=0,04) foi rejeitado atomicamente com 6 passos aceitos (0,24 min; voxel 769), enquanto cB_q043 (dt=0,02) completou 2 min sem rejeicao. A 0,43 L/min a extracao de O2 e ~20%, o O2 dissolvido do tecido cai a ~0,065-0,07 mM
e o estoque do voxel (~5e-4 pmol) fica abaixo do consumo de 3 hospedeiros em 0,04 min (5,5e-4 pmol): limite de separacao de operadores, NAO escassez fisiologica. Portanto a P6 NAO foi testada por esse caso (nao e confirmacao nem refutacao).
Caso de reposicao: cE2_q043_long (Q=0,43, dt=0,02, beta=107, 600 passos = 12 min), a rodar apos os demais casos terminarem (evitar contencao de CPU). A P6 vale para ele sem alteracao: rejeicao por escassez entre 6 e 20 min; se nao houver rejeicao em 12 min, a P6 fica sem desfecho (so cota inferior).
Risco declarado: cE_q086_long (dt=0,04, Q=0,86) opera perto do mesmo limite (razao consumo/estoque 0,98 no segmento central com 3 hospedeiros); uma rejeicao nele deve ser lida como artefato de passo, nao escassez.
Proposta registrada (NAO executada): acoplar ao LP o estoque local de sangue do voxel (parcela do conteudo do no, com Hb) removeria o limite de passo; exige previsoes proprias antes de rodar.

## ADENDO 4: acoplamento do estoque local de sangue (LBP), registrado ANTES de rodar (binario native_c2, derivado de native_c por make_native_c2.py)
Motivo: o limite de passo (dt <~ 0,04 min) vem de o LP enxergar so o O2 dissolvido do tecido do voxel; fisicamente o voxel tem sangue (com Hb) a seu lado. Mudanca de MODELO declarada (hipotese exploratoria, nao evidencia): com local_pool=1 o LP ve o estoque do tecido + a parcela do voxel no conteudo do seu no sanguineo
(conteudo do no / numero de voxels do no); a captacao alem do estoque do tecido e debitada do no (conserva massa; o tecido nunca fica negativo; so um reajuste de arredondamento < 1e-12 pmol e permitido e gravado em lbp.json). Sem local_pool o comportamento e o de native_c.
Casos (todos beta=107,45, local_pool=1; passo e horizonte): cL1 (Q=10, dt=0,1, 20 passos), cL2 (Q=10, dt=1,0, 10 passos), cL5 (Q=0,86, dt=0,02, 100 passos; comparar com cB_q086 sem LBP), cL6 (Q=0,86, dt=0,2, 25 passos = 5 min), cL3 (Q=0,86, dt=0,1, 300 passos = 30 min),
cL4 (Q=0,43, dt=0,1, 300 passos = 30 min), cL7 (Q=0,86, dt=0,05, 100 passos = 5 min); o refinamento a 5 min compara cL6, cL3 (passo 50) e cL7.
Previsoes (escritas antes de rodar):
PL1 (alta): cL1 e cL2 reproduzem o cenario estatico (B0b) a <= 1e-6 por agente em todos os passos, balancos <= 1e-9, sem rejeicao (o caso equivalente sem LBP a dt=0,1 foi rejeitado).
PL2 (alta): cL5 e cB_q086 aplicam fluxos identicos (diferenca 0,0) em 2 min; os estoques armazenados diferem < 2%.
PL3 (media): em Q=0,86 (cL3) a captacao de glicerol cai abaixo do teto (1,38 pmol/min) em ate 10 min e fica <= 1,19 pmol/min (entrada) no fim; o lactato sobe acima de 8,02 pmol/min; sem rejeicao em 30 min.
PL4 (media): em Q=0,43 (cL4) ha rejeicao atomica por escassez entre 5 e 30 min (a oferta de carbono gliconeogenico limitada: lactato 1,10 + glicerol <=0,10 + alanina 0,12 unidades dao <~0,66 de glicose gliconeogenica contra 0,86 necessaria, o que o estoque so adia).
PL5 (media): no refinamento a 5 min (dt 0,2; 0,1; 0,05) as captacoes totais e os estoques mudam < 1% por divisao do passo em todos os campos com estoque > 1e-12 pmol, exceto talvez glicerol.
PL6 (alta): lbp.json mostra reajuste de arredondamento < 1e-12 pmol.
Se alguma falhar, reportar. Nenhum parametro sera ajustado.

## ADENDO 5 (2026-11-20, apos o reinicio do contêiner)
O contêiner reiniciou por volta de 05:26 e interrompeu cE2_q043_long (36 de 600 passos, sem resultado), cL5 e cL6/cL7/cL3/cL4 (nao concluidos); nada disso gerou conclusao. cE2 NAO foi relancado (cL4, com LBP, cobre Q=0,43; a P6 sem LBP fica sem desfecho alem da cota inferior de 0,8 min do cB_q043 em 2 min).
Resultados ja obtidos antes do reinicio (nao afetados): cL1 (LBP, Q=10, dt=0,1, 20 passos): desvio 0,0 do estatico, balanco 1,5e-12, reajuste de arredondamento 0 -> PL1 (parte dt=0,1) e PL6 CONFIRMADAS.
cL2 (LBP, Q=10, dt=1,0): REJEITADO atomicamente no primeiro passo (voxel 67) -> PL1 FALHOU para dt=1,0: o LBP estende o passo util de ~0,04 para ~0,1 min mas nao para 1 min (o conteudo do no, ~0,35 pmol de O2 para 4-205 voxels, e menor que o consumo de um minuto dos voxels do no; a reposicao por fluxo dentro do passo nao e vista pelo LP).
Relancamento: cadeia cL6 -> cL5 -> cL7 e, em paralelo, cL3 e cL4 (run_cases_v5b.sh).

Nota (2026-11-20): cL6 (LBP, Q=0,86, dt=0,2) foi REJEITADO atomicamente (rc=42): o LBP nao sustenta dt=0,2; o refinamento a 5 min com LBP fica restrito a dt=0,1 e 0,05 (cL3 passo 50 e cL7). Autorizacao do pesquisador para refazer o que se perdeu no reinicio: cE2_q043_long
(Q=0,43, dt=0,02, sem LBP, 600 passos) foi enfileirado para depois do lote v5b (limite de memoria ~3,8 GB com 3 execucoes simultaneas).
