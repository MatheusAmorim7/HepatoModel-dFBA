# Checagem de ordem de grandeza contra dados humanos independentes (2026-11-22)
Exploratória; NÃO é validação nem calibração. Valores humanos vêm de resumos de busca na web (não li os artigos completos); conferir os originais antes de citar. Modelo: Z0 (controle), 70 kg, 1 µmol/min = 0,0031194 unidades (anchors.py), 30 min estacionário; E3_s3_k10 e E3_s3_k15 entre parênteses.

| Observável | Modelo | Humano (fonte) | Leitura |
|---|---|---|---|
| Produção de glicose | 9,8 µmol/kg/min (9,7; 11,2) | ENTRADA do modelo (âncora), não independente | — |
| Fração gliconeogênica | 40% (piso imposto) | 23–42% após 14 h; 36% (0,79 de 2,19 mg/kg/min) após jejum noturno (Landau, JCI; PMC508308) | Dentro, mas no limite superior; é a nossa escolha de piso, não resultado |
| Glicose gliconeogênica | ≈3,9 µmol/kg/min | ≈4,4 (36% de 2,19 mg/kg/min) | Compatível |
| Glicose derivada de lactato | ≈3,2 µmol/kg/min (lactato 6,3 captado) | 1,97 µmol/kg/min sistêmico (fígado + rim; Meyer 2002) | Modelo ~1,6× acima do total sistêmico; mistura de precursores muito rica em lactato (78%) |
| Consumo de O₂ do fígado | **16,7 mL/min** (16,5; 19,7) | esplâncnico 40,7 mL/min/m² (Kameda) ≈ 73 mL/min a 1,8 m² (inclui vísceras portais; fígado sozinho é menor) | **Modelo ~2,5–4× abaixo**: maior discrepância |
| Captação de palmitato | 30 µmol/min | oxidação hepática de AG 0,08–0,16 µmol/min/mL (Iozzo 2010; unidades e volume hepático não reverificados) | Provável subestimativa (não conclusivo) |

## Leitura
1. **Compatível:** fração gliconeogênica e glicose gliconeogênica total (por construção parcial: piso 40%).
2. **Alta:** o lactato como precursor (~1,6× o dado sistêmico) e o teto de glicerol (0,86 µmol/kg/min) foram definidos por âncoras/tetos, não medidos aqui.
3. **Discrepância relevante: O₂.** O fígado do modelo consome 3 a 4× menos O₂ que o esplâncnico medido. Consequência direta: a queda de O₂ no sinusoide é de apenas ~4–10% no modelo, enquanto in vivo o gradiente periportal→pericentral é grande. O ATPM=5 (único sorvedouro de ATP não metabólico) provavelmente subestima o gasto energético (síntese de proteína, bombas iônicas, VLDL), o que também pode explicar o palmitato baixo.
4. **Implicação:** o gradiente de O₂ é o principal motor biológico da zonação e da hipóxia local; sem ele o modelo não discrimina zonas por oxigênio. Isto limita a interpretação dos mapas de O₂ e do tumor.
## Próximo passo sugerido
Calibrar o ATPM (ou acrescentar um sorvedouro de ATP) contra o consumo hepático de O₂ como alvo EXPLÍCITO, com previsões registradas antes (viabilidade do cenário, lactato, palmitato, gradiente de O₂), mantendo a versão atual como controle. Não feito: exige escopo novo e vira calibração declarada.
