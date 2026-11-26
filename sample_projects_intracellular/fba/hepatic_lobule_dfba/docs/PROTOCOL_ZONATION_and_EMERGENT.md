# Protocolo Z: zonacao enzimatica no PhysiCell nativo (registrado ANTES de rodar qualquer caso)

Escopo (autorizacao: "Isso, vamos concluir essa etapa 1..."): copia parametrizada do nativo da Fase C (native_c2 -> native_z), sem tocar arquivos selados nem a pasta da Fase C, com tarefa de glicose e tetos enzimaticos POR POSICAO periportal->central.
Fora de escopo: tumor, ciclo/morte/movimento, zonacao de ureia/GLUL/GK/PK, mudanca do modelo raiz/Git, parametros fisiologicos adotados, promocao formal de estagio. Exploratorio; nao e calibracao nem validacao.

## Desenho (hipoteses, nao medidas)
- Eixo x=(segmento+0,5)/8 do voxel de cada um dos 3000 hospedeiros (cells.csv; segmentos 0=periportal ... 7=central; hospedeiros por segmento 700,611,508,430,341,221,138,51).
- Abundancia PP-alta (PEPCK r0123, FBP1 r0487, G6Pase r0396): E(x)=1+(1/s-1)x, normalizada para media 1 sobre os HOSPEDEIROS (m). s=1 controle; s=3 (unico numero com fonte: PK 3x no rato, Berndt 2018) e s=10 (sensibilidade).
- TAREFA por hospedeiro: exportacao de glicose 2,14*m (lb=ub) e teto de GLYCOGEN_SUPPLY *m (preserva >=40% gliconeogenese); media sobre hospedeiros = 2,14 => total do orgao igual ao controle.
- TETOS: ub = kappa*m*J_ref (kappa=1,5; J_ref = fluxo estatico por hospedeiro: PEPCK 1,482, FBP1 0,856, G6Pase 2,0233 do B0b).
- Por que tarefa por posicao: o S9 mostrou que sem tarefa no nivel do orgao o problema e degenerado; no nativo cada hospedeiro resolve o LP local, entao a divisao do trabalho entre posicoes tem de entrar como tarefa. A gradiente de fluxo gliconeogenico e portanto PRESCRITA (nao emergente); o que e emergente: gradientes de substrato, extracao, viabilidade local, estado estacionario.
- Nao zonados (e por que): ureia/CPS1 (S9: PP:PC nao robusto; exigiria escalar tetos de N), GLUL (fora do pacote v1), GK/PK (fluxo 0 no cenario pos-absortivo).
## Bracos (Q=0,86 L/min, dt=0,1 min, estoque local de sangue, 300 passos=30 min; configuracao do cL3)
Z0_uniform (todos os fatores 1; controle) | Z3_s3_k15 | Z10_s10_k15 | Z3c_caps_only (s=3 com tetos, SEM tarefa por posicao: controle negativo)
## Previsoes (escritas antes de rodar)
P1 (alta): Z0 reproduz o cL3 (fluxos aplicados por hospedeiro e passo) a <=1e-6 e fecha balancos <=1e-9 pmol.
P2 (alta): Z3c e rejeitado atomicamente (LP local inviavel) no passo 0 em voxel de segmento >=4, porque o teto de PEPCK/FBP1/G6Pase no pericentral (kappa*m*J_ref, m~0,5) fica abaixo da tarefa uniforme.
P3 (alta): Z3 e Z10 viaveis por 30 min; totais do orgao (captacao de lactato, glicerol, O2; exportacao de glicose) diferem de Z0 <=3% no estado quase estacionario (t>=10 min); balancos <=1e-9.
P4 (alta): fluxo de PEPCK por segmento em Z3 decresce monotonicamente PP->PC com razao seg0/seg7 dentro de 10% de 1,2353/0,4834=2,56 (confere a prescricao).
P5 (media): em Z3 a queda PP->PC da concentracao de lactato no sangue e MAIOR que em Z0 (Z0: -21% no cL3); Z10 > Z3.
P6 (media): a extracao de glicerol (>=0,9 em Z0) nao cai mais que 0,03 em Z3; a captacao de glicerol se concentra no PP (fracao do total nos segmentos 0-2 maior em Z3 que em Z0).
P7 (media): estado quase estacionario em <=6 min; deriva <1% entre 20 e 30 min em todos os bracos viaveis.
Se falhar, reportar; sem ajustar s, kappa ou pesos.

# ADENDO E (2026-11-22, registrado ANTES de rodar): fluxo gliconeogenico EMERGENTE
Pedido do pesquisador: "podemos tentar emergir esse fluxo". Mudanca: sem tarefa de glicose por hospedeiro. Cada hospedeiro MAXIMIZA a exportacao de glicose (recompensa lambda=200 por unidade >> custo L1 marginal; exportacao livre em [0,10]) sujeito a capacidade, suprimento (estoques/sangue), ureia 0,19 e ATPM 5.
Entrada passa a ser a ABUNDANCIA ENZIMATICA (tetos de PEPCK/FBP1/G6Pase = kappa*m*J_ref; teto de glicogenolise *m; hipotese linear s=3), e o fluxo e previsao. Binario native_ze (recompensa so no custo do LP; recompensa 0 => identico ao native_z).
Bracos: E0_uniform_k10 (s=1,kappa=1) | E3_s3_k10 | E0_uniform_k15 (s=1,kappa=1,5) | E3_s3_k15. Q=0,86, dt=0,1, estoque local, 30 min. kappa=1 => capacidade = referencia pos-absortiva (exportacao max 2,14); kappa=1,5 => capacidade excedente, suprimento decide.
Previsoes (antes de rodar):
Q1 (alta): E0_k10 reproduz a referencia: exportacao de glicose por hospedeiro 2,14 (+-0,02), totais do orgao <=2% de Z0, viavel 30 min.
Q2 (media-alta): E3_k10 viavel, exportacao media 2,14 (+-3%) e PEPCK/FBP1 no TETO em todos os segmentos (fluxo/teto >=0,99): sem suprimento limitante, o gradiente de fluxo = gradiente de capacidade (sem zonacao nova emergente).
Q3 (media): E3_k15: exportacao media do orgao >2,14 (+10 a +40%), extracao de lactato >0,8 e segmentos >=5 abaixo do teto (fluxo/teto <0,9): gradiente de fluxo mais ingreme que o de capacidade (razao PEPCK seg0/seg7 > 2,56).
Q4 (media): E0_k15 (enzimas uniformes): ainda emerge gradiente PP>PC so por suprimento (PEPCK seg0/seg7 >=1,2) e exportacao media >2,14.
Q5 (media-baixa): estado quase estacionario em <=6 min e deriva <1% nos bracos viaveis (hospedeiros gananciosos podem oscilar; relatar se nao).
Q6 (alta): balancos <=1e-9, sem rejeicao.
