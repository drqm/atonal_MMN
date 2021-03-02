(start-idyom)
(idyom-db:describe-database :verbose t)
(idyom-db:describe-dataset 23 :verbose t)

;(idyom:idyom 26 '(cpitch) '(cpitch) :models :stm :k 1 :detail 3 :output-path "/media/sf_Documents/projects/atonal_mumufe_beh/")
;(idyom:idyom 11 '(cpitch) '(cpitch) :models :stm :k 1 :detail 3 :output-path "/media/sf_Documents/projects/atonal_mumufe_beh/")

(idyom-db:import-data :mid "/media/sf_Documents/projects/atonal_mumufe_beh/midi_corpus_pclass/24/" "Roger Dean's algorithmic serial melodies - pclass" 34)

(idyom-db:import-data :mid "/media/sf_Documents/projects/atonal_mumufe_beh/midi_corpus_pclass/25/" "Hindemith, Varese - Pieces for solo flute - pclass" 35)

(idyom-db:import-data :mid "/media/sf_Documents/projects/atonal_mumufe_beh/non-tonal_corpora/WebernPianoVariations/" "Webern Piano Variations" 26)

(idyom-db:import-data :mid "/media/sf_Documents/projects/atonal_mumufe_beh/midi_corpus_pclass/10/" "Songs and ballads from Nova Scotia, Canada, transposed to c and pclass" 30)

(idyom-db:import-data :mid "/media/sf_Documents/projects/atonal_mumufe_beh/midi_corpus_pclass/11/" "Chorale melodies harmonised by J.S. Bach, transposed to c and pclass" 31)

(idyom-db:import-data :mid "/media/sf_Documents/projects/atonal_mumufe_beh/midi_corpus_pclass/12/" "German folk songs from the Essen Folk Song Collection: fink, transposed to c and pclass" 32)

(idyom-db:import-data :mid "/media/sf_Documents/projects/atonal_mumufe_beh/midi_mels_pclass/tonal/" "tonal melodies pclass (atonal mumufe project)" 28)

(idyom-db:import-data :mid "/media/sf_Documents/projects/atonal_mumufe_beh/midi_mels_pclass/atonal/" "atonal melodies pclass (atonal mumufe project)" 29)

;(idyom-db:import-data :mid "/media/sf_Documents/projects/atonal_mumufe_beh/midi_recordings/" "atonal recordings (atonal mumufe project)" 23)

;(idyom:idyom 8 '(cpitch) '((cpitch-class cpcint)) :models :both+ :k 1 :pretraining-ids '(20 21 22) :detail 3 :output-path  "/media/sf_Documents/projects/atonal_mumufe_beh/IDyOM_models/")

;(idyom:idyom 9 '(cpitch-class) '((cpitch-class cpcint)) :models :stm :k 1 :detail 3 :output-path  "/media/sf_Documents/projects/atonal_mumufe_beh/IDyOM_models/")


(idyom:idyom 28 '(cpitch) '(cpitch) :models :both+ :k 1 :pretraining-ids '(30 31 32) :detail 3 :output-path  "/media/sf_Documents/projects/atonal_mumufe_beh/IDyOM_models_pclass/")

(idyom:idyom 29 '(cpitch) '(cpitch) :models :both+ :k 1 :pretraining-ids '(30 31 32) :detail 3 :output-path  "/media/sf_Documents/projects/atonal_mumufe_beh/IDyOM_models_pclass/")

(idyom:idyom 28 '(cpitch) '(cpitch) :models :both+ :k 1 :pretraining-ids '(34 35) :detail 3 :output-path  "/media/sf_Documents/projects/atonal_mumufe_beh/IDyOM_models_pclass/")

(idyom:idyom 29 '(cpitch) '(cpitch) :models :both+ :k 1 :pretraining-ids '(34 35) :detail 3 :output-path  "/media/sf_Documents/projects/atonal_mumufe_beh/IDyOM_models_pclass/")

;(idyom:idyom 28 '(cpitch) '((cpitch cpint)) :models :both+ :k 1 :pretraining-ids '(23) :detail 3 :output-path  "/media/sf_Documents/projects/atonal_mumufe_beh/IDyOM_models_pclass/")

;(idyom:idyom 29 '(cpitch) '((cpitch cpint)) :models :both+ :k 1 :pretraining-ids '(23) :detail 3 :output-path  "/media/sf_Documents/projects/atonal_mumufe_beh/IDyOM_models_pclass/")
'
(idyom:idyom 28 '(cpitch) '(cpitch) :models :stm :k 1 :detail 3 :output-path  "/media/sf_Documents/projects/atonal_mumufe_beh/IDyOM_models_pclass/")

(idyom:idyom 29 '(cpitch) '(cpitch) :models :stm :k 1 :detail 3 :output-path  "/media/sf_Documents/projects/atonal_mumufe_beh/IDyOM_models_pclass/")
