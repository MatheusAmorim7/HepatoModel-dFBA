# O que o executavel nativo atual fixa e o que a Fase B precisa parametrizar

Fonte: notebooks/lobule_metabolic_cadence_2026_11_17/native.cpp, cache.h, workflow.py (selados; NAO modificar). A Fase B usa uma COPIA em pasta nova.

## Fixo no codigo atual (contagens por grep; numeros dentro de native.cpp)
- 179 reacoes e 175 metabolitos: `need(...reactions.size()==179 && ...metabolites.size()==175, "Individual_model_shape")`; laços `for j<179`; arrays 3000x179 (bounds.f64, fluxes).
  Pacote v1: 225 reacoes, 217 metabolitos.
- 20 campos de recurso (estoques compartilhados por voxel): `fields.size()==20`, `initial_stocks` 2500x20 no XML (`finite_resource_control`), laços `f<20`, LP conjunto por voxel com
  linhas de estoque (`stock_start+f`) e blocos `c*(2*n+20)`. O pacote exige 25 campos: + palmitato (EX_HC00226_s), arginina, ornitina, citrulina e glicerol (EX_HC00114_s);
  GLYCOGEN_SUPPLY NAO e campo (e limite por agente).
- 3000 agentes e 2500 voxels, 1884 voxels ocupados (`need(occupied==1884)`), contagens 537000/525000 (3000x179, 3000x175), `fields[0]==oxigenio`, conversao `*.12*dt/60`, volume do voxel `.008`.
- Metabolismo dentro de PhysiCell: uma rede SBML individual por agente (libRoadRunner-like via dFBA addon; `CoordinatedHost`), selecao L1 global por voxel (Clp, tolerancias 1e-10), cache de LP por
  voxel com membros e limites fixos, rejeicao sem reparo, escrita comprimida. Passos testados: 0,0005 min; horizonte 201 passos (6,03 s simulados).
## O que o cenario de referencia pos-absortivo exige do nativo
- Todas as restricoes do cenario sao LIMITES por reacao (scenario_reference_f1.csv / f2.csv), inclusive a particao de glicogenolise (GLYCOGEN_SUPPLY ub = 0,6 G / 11). Nenhuma linha extra de LP.
- Igualdades por agente: ureia [0,19; 0,19] e glicose [2,14; 2,14]. A selecao L1 deve aceitar limites inferior=superior.
- Os cinco tetos que variam por agente (r0353, r0396, r0123, r0487, r0054) vem de source.npz; a coluna de limites por agente tem que cobrir as 225 reacoes.
- Estoques: para o teste de paridade, estoques amplos (nao limitantes) reproduzem a solucao estatica; so na Fase C a oferta e o transporte vascular os tornam limitantes.
## Teste de paridade proposto (B0)
1 passo, 3000 agentes, mesmos limites do cenario por agente, estoques nao limitantes: fluxos do LP conjunto nativo (Clp) vs LP em Python do mesmo pacote; criterio |diff| <= 1e-6 nos 25 campos de troca e
no objetivo L1 por voxel; Sv, limites e balanco C/N conferidos como nos controles anteriores. Sem tumor, sem zonacao, sem ciclo/morte/movimento.
## Governanca
Trabalho nativo novo exige aprovacao autonoma de escopo exato (artefato protocol.json com hashes, como a de 2026-11-17); nenhuma aprovacao antiga cobre estes arquivos. Nenhuma sessao formal
e criada. O agente PhysiCell (configurador) so entra quando houver configuracao de dominio/substratos/regras a gerar (Fase C/D).
