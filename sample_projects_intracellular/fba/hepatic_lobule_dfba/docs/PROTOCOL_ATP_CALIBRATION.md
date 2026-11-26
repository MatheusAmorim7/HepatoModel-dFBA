# Protocolo ATP: calibração declarada do ATP de manutenção contra o consumo hepático de O2 (registrado ANTES de rodar no nativo)

Motivação: CHECAGEM_EXTERNA_PT.md (lobule_zonation_native_2026_11_21): O2 hepático do modelo 16,7 mL/min vs alvo humano; gradiente de O2 no sinusoide ~4–10%. Isto é CALIBRAÇÃO DECLARADA de um único parâmetro (CULTURE_ATPM) contra UM alvo; o resto é validação não usada no ajuste. Exploratório.
## Alvo e parâmetro
- Alvo de ajuste: consumo hepático de O2 = 42,9 mL/min = Elia/Wang 2010 (200 kcal/kg/dia, taxa específica do fígado; resumo de busca, original não lido) × 1,5 kg de fígado (âncora do modelo: 107e6 hepatócitos/g × 21,4 g/kg × 70 kg) ÷ 4,86 kcal/L O2 (equivalente energético assumido, RQ~0,8). Faixa de literatura mais alta (6 mL/100 g/min × 1500 g = 90 mL/min; Lutz 1975 via resumo secundário) NÃO é alvo; é analisada como limite estrutural.
- Parâmetro ajustado: CULTURE_ATPM (lb=ub) de 5 para 21,6 (interpolação estática L1 em scan_atpm.py: 42,94 mL/min; SCAN_ATPM.json). Ajuste feito ANTES de rodar no nativo, sem olhar outros observáveis.
- Mudança declarada auxiliar: teto de captação de O2 por hospedeiro de −6 para −20 (o −6 era uma constante antiga; truncaria hospedeiros periportais zonados; O2 real é limitado pelo transporte).
- Mudança declarada de especificação (não ajuste): em braços "pal045" o campo de ácidos graxos (modelado só como palmitato C16) entra com o equivalente em carbono do NEFA medido (CHS: palmítico 124, esteárico 60, oleico 150, linoleico 79 µmol/L; ×18/16 para C18 => 449 µM), em vez de só 124 µM do palmitato.
## Bracos (Q=0,86 L/min, dt=0,1, estoque local, 30 min; native_z)
A1_uniform_pal124 (ATPM 21,6, sem zonação, palmitato 0,124) | A2_uniform_pal449 (idem, 0,449 mM) | A3_s3_pal449 (idem A2 com zonação prescrita s=3 do Z3)
## Previsões (escritas antes de rodar)
R1 (alta): estático, ATPM=21,6 dá 42,9±1 mL/min (já visto no scan; usar como verificação) e o O2 do órgão no nativo (A2) fica em 42,9±5%.
R2 (média-alta): A1 (palmitato só 0,124 mM) é rejeitado atomicamente (LP local inviável) antes de 10 min por falta de palmitato (a extração média exigida ≈0,75 seria satisfeita só no PP).
R3 (média): A2 viável por 30 min; a queda do conteúdo de O2 no sangue do segmento 0 ao 7 passa de ~4% (ATPM=5) para ≥20%.
R4 (média): captação hepática total de palmitato em A2 = 81±20 µmol/min, abaixo da faixa derivada de Iozzo (120–240 µmol/min; unidades não reverificadas): previsão de discrepância residual (validação não ajustada).
R5 (média-baixa): totais de lactato, glicerol e glicose do órgão em A2 diferem ≤5% do controle ATPM=5 (tarefas fixam a gliconeogênese).
R6 (média): A3: captação de O2 por hospedeiro periportal/central ≥1,4 e nenhum hospedeiro rejeitado por O2 (teto −20 não binding).
R7 (alta): balanços por campo ≤1e-9 pmol, sem reparos.
Validação NÃO ajustada (reportar sem ajustar nada): palmitato vs Iozzo, extração hepática de O2 (~0,37, Lutz via resumo secundário), perfil PP→PC do conteúdo de O2, lactato.
Se falhar, reportar; sem ajustar ATPM de novo.

# ADENDO A (2026-11-22, registrado ANTES de rodar): binário com contrato de ATPM parametrizado
Fato (observado na primeira tentativa, 3 braços rc=1): o nativo tem a verificação fixa `ATPM contract` (CULTURE_ATPM==5) em src/zonation_2d_native_groups/native.cpp:71; o escopo registrado dizia "binário native_z inalterado", o que era inviável. Correção: cópia própria da árvore de fontes na pasta ATP com a verificação lendo ATPM_EXPECTED do ambiente (padrão 5 => comportamento idêntico). Binário `native_zatp` (make_native_zatp.py); native_z e a pasta da zonação intactos. Os braços rodam com ATPM_EXPECTED=21.6. Nenhuma previsão muda.

# ADENDO B (2026-11-22, escrito DEPOIS de ver a falha dos 3 braços a dt=0,1; resultados preservados em runs/failed_dt01_step0)
Fato: com ATPM=21,6 os três braços foram rejeitados atomicamente no passo 0 (local_LP_infeasible, voxel 69, 3 hospedeiros): o consumo de O2 de um passo de 0,1 min supera o O2 disponível (tecido + parcela do sangue do nó), o mesmo limite de separação de operadores visto na Fase C (consumo/estoque), agora 2,6x pior. Não é resultado científico.
Mudança declarada (pos-hoc): dt=0,02 min (dt_mechanics 0,1 = 5 passos; dt_phenotype 2,0) e horizonte 6 min (300 passos), em vez de dt=0,1 e 30 min; os estados estacionários anteriores se estabeleceram em ≤4 min. R7 inalterada; R3/R5/R6 avaliadas ao fim de 6 min (nao 30). Se ainda houver rejeição, relatar.
