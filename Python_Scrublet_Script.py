samples = ['HPAP-019','HPAP-020','HPAP-021','HPAP-022','HPAP-023','HPAP-024','HPAP-026','HPAP-028','HPAP-029','HPAP-032','HPAP-034','HPAP-035','HPAP-036','HPAP-037','HPAP-038','HPAP-039','HPAP-040','HPAP-042','HPAP-043','HPAP-044','HPAP-045','HPAP-047','HPAP-049','HPAP-050','HPAP-051','HPAP-052','HPAP-053','HPAP-054','HPAP-055','HPAP-056','HPAP-057','HPAP-058','HPAP-059','HPAP-061','HPAP-063','HPAP-064','HPAP-065','HPAP-070','HPAP-071','HPAP-072','HPAP-074','HPAP-075','HPAP-077','HPAP-079','HPAP-080','HPAP-081','HPAP-082','HPAP-083','HPAP-084','HPAP-085','HPAP-087','HPAP-088','HPAP-090','HPAP-091','HPAP-092','HPAP-099','HPAP-100','HPAP-101','HPAP-103','HPAP-104','HPAP-105','HPAP-106','HPAP-107','HPAP-108','HPAP-109']
for sample in samples:
    #Read in sample counts and gene lists
    input_dir = '~/hpap/Scrublet/' + sample
    counts_matrix = scipy.io.mmread(input_dir + '_matrix_roundToInt.mtx').T.tocsc()
    genes = np.array(scr.load_genes(input_dir + '_genes_roundToInt.tsv', delimiter='\t', column=0))
    
    #Run Scrublet with default thresholds
    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                              min_cells=3, 
                                                              min_gene_variability_pctl=85, 
                                                              n_prin_comps=30)
    
    #Pull out predicted doublets from Scrublet
    predicted_doublets = scrub.call_doublets() 
    
    #Save results
    barcodes = input_dir + '_barcodes_roundToInt.tsv'
    with open(barcodes, 'r') as f:
        with open('~/hpap/Scrublet/scrublet_predicted_doublets_cutoff{}.txt'.format(scrub.threshold_), 'w') as o:
            counter=0
            for line in f:
                o.write('\t'.join((line.strip(), str(predicted_doublets[counter]), str(doublet_scores[counter]), '\n')))
                counter+=1
    
    #Delete the Scrublet object to continue loop
    del scrub
