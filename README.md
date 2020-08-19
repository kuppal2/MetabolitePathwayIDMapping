# Metabolite Pathway ID Mapping

By: Karan Uppal

***METHODS / WORKFLOW***

- R packages used: webchem, KEGGREST, RJSONLITE
- Online sources: PUG (PubChem) REST, KEGG REST, Chemical Translation Service
- Convert PubChem CID to KEGG
   + Convert PubChem CID to PubChem SID using PUG REST
   + Convert PubChem SID to KEGG ID using keggConv function
   + If no matches, then
      - Get InChiI keys using cts_convert function in R package webchem
      - Get ChEBI ID using get_chebiid function using InChi keys
      - Convert ChEBI ID to KEGG ID using keggConv function
- Convert PubChem CID to ChEBI ID
   + Get InChiI keys using cts_convert function in R package webchem
   + Get ChEBI ID using get_chebiid function using InChi keys
- Download KEGG pathway molecular atlas 
   + Use the keggList function to get pathway IDs for input species code (e.g. has”)
   + For each pathway ID, use the keggGet function to get genes and compounds associated with that pathway
   + Map the KEGG gene IDs to NCBI Gene IDs and NCBI Protein IDs using the keggConv function
- Download Reactome pathway molecular atlas 
   + Download and merge the following files:
      - https://reactome.org/download/current/NCBI2Reactome_PE_All_Levels.txt
      - https://reactome.org/download/current/ChEBI2Reactome_PE_All_Levels.txt
      - https://reactome.org/download/current/UniProt2Reactome_PE_All_Levels.txt
      - Filter by species code (e.g. ”Homo sapiens”)
