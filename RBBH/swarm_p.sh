blastp -task blastp -query scap.faa -db shom_p -outfmt 6 -out CvH.blastp
blastp -task blastp -query shom.faa -db scap_p -outfmt 6 -out HvC.blastp
blastp -task blastp -query shom.faa -db sepi_p -outfmt 6 -out HvE.blastp
blastp -task blastp -query sepi.faa -db shom_p -outfmt 6 -out EvH.blastp
blastp -task blastp -query sepi.faa -db scap_p -outfmt 6 -out EvC.blastp
blastp -task blastp -query scap.faa -db sepi_p -outfmt 6 -out CvE.blastp
