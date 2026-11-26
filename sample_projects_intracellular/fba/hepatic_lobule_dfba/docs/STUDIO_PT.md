# Como ver as simulações no PhysiCell Studio

Cada pasta em `studio_output/` é uma rodada convertida para a saída padrão do PhysiCell (XML + .mat). Validada com o leitor do próprio Studio (`pyMCDS` / `pyMCDS_cells` em PhysiCell-Studio-dfba/bin): 25 substratos e 3000 células carregam sem erro. NÃO testei a interface gráfica do Studio (não há tela aqui).
Pastas: Z0_uniform, Z3_s3_k15 (zonação prescrita), E3_s3_k10 e E3_s3_k15 (fluxo emergente), A2_uniform_pal449 e A3_s3_pal449 (ATPM calibrado, O₂ hepático 43 mL/min).
Snapshots 0–5: passos 0, 15, 30, 75, 150, 299 (tempo = (passo+1)×dt; dt=0,1 min nas pastas Z/E e 0,02 min nas A).
## Abrir
No Studio (`python bin/studio.py -c config/PhysiCell_settings.xml`, o config é só para abrir): aba **Plot** → campo do diretório de saída → apontar para uma das pastas acima e escolher `output0000000N.xml`. O Studio não tem opção de linha de comando para a pasta de saída (conferido em studio.py); o diretório se escolhe na aba. Se a aba Plot reclamar de substratos que não estão no config, use a aba Plot só com o diretório de saída ou me avise que gero um config com os 25 substratos.
- Substratos (menu da aba Plot): os 25 campos em mM (concentração no tecido após o passo). Bom para: lactate, glycerol, oxygen, ammonia, palmitate.
- Células: cor por `cell_type` (zone_1/2/3), `segment_PP0_to_PC7` (0 periportal … 7 central) ou pelas variáveis de fluxo: PEPCK_r0123, FBP1_r0487, G6Pase_r0396, glycogenolysis, glucose_export, lactate_flux, glycerol_flux, O2_flux, urea_export, palmitate_flux, alanine_flux, NH3_flux (unidades de fluxo do modelo; negativo = captação nas trocas).
## Cuidados
- É só visualização do que o binário nativo calculou; não é resultado novo nem validação.
- Valores por célula dentro de um mesmo voxel são uma divisão degenerada do LP (o total do voxel é determinado, a divisão entre hospedeiros não): leia agregados por voxel ou segmento.
- Sem difusão, ciclo, morte ou movimento: as células são fixas e o "transporte" é o sangue serial, não a difusão do PhysiCell; o Studio não reproduz a simulação, só mostra a saída.
- Para reexportar: `python export_physicell.py` (usa as rodadas em notebooks/lobule_zonation_native_2026_11_21/runs e lobule_atp_calibration_2026_11_22/runs).
