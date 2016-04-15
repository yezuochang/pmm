function ncss_export(G,W,subcktname,outputfile,opts)
export_spectre_model(G,ss_transpose(W), subcktname, outputfile,'w');