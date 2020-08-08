library(webchem)
library(KEGGREST)
library(plyr)
library(RJSONIO)


#input: PubChem CIDs
#output: data frame with PubChem CIDs mapped to KEGG IDs and PubChem SIDs (Substance IDs)
pubchem_to_kegg<-function(cids){

	cids<-na.omit(cids)

	if(length(cids)>10){
		
		res_temp<-lapply(seq(1,length(cids),10),function(j){

			cur_range<-j:(j+9)
			
			res<-pubchem_to_kegg_child(cids[cur_range])
            res<-unique(res)
            return(res)
			})
        
		res_final<-ldply(res_temp,rbind)
	}else{
	
		res_final<-pubchem_to_kegg_child(cids)

	}
    colnames(res_final)<-c("pubchem_cid","pubchem_sid","KEGGID")
	return(res_final)

}

pubchem_to_kegg_child<-function(cids){


	cids<-na.omit(cids)
	

	if(length(cids)>0){
	
	#PUG REST JSON	
	pug_json<-paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",paste(cids,collapse=","),"/sids/JSON",sep="")
	
	d1=fromJSON(pug_json)

	if(length(d1)>0){

	sid_list<-lapply(1:length(d1[[1]]$Information),function(j){
        

		if(length(d1[[1]]$Information[[j]])>1){
            Sys.sleep(0.1)
            kegg_ids<-keggConv(target="compound", source=paste("pubchem:",tolower(d1[[1]]$Information[[j]]$SID),sep=""),querySize = 100)
            
            sid_list<-paste(d1[[1]]$Information[[j]]$SID,collapse=";")
            
            if(length(kegg_ids)<1){
                
                    chebi_ids<-pubchem_to_chebi(d1[[1]]$Information[[j]]$CID)
                    chebi_ids<-tolower(chebi_ids$chebiid)
                    kegg_ids<-keggConv(target="compound", source=chebi_ids,querySize = 100)
                    
                    
                     resmat<-as.data.frame(cbind(names(kegg_ids),kegg_ids))
                           resmat<-unique(resmat)
                           if(length(kegg_ids)>0){
                               resmat$V1<-(paste(d1[[1]]$Information[[j]]$SID,collapse=";"))
                               
                               colnames(resmat)<-c("pubchem_sid","KEGGID")
                               resmat$pubchem_sid<-gsub(resmat$pubchem_sid,pattern="pubchem:",replacement="")
                               resmat$KEGGID<-gsub(resmat$KEGGID,pattern="cpd:",replacement="")
                               rownames(resmat)<-NULL
                               
                               
                               res_final<-cbind(d1[[1]]$Information[[j]]$CID,resmat)
                               colnames(res_final)<-c("pubchem_cid","pubchem_sid","KEGGID")
                           }else{
                               
                                   res_final<-cbind(d1[[1]]$Information[[j]]$CID,(paste(d1[[1]]$Information[[j]]$SID,collapse=";")),NA)
                                   colnames(res_final)<-c("pubchem_cid","pubchem_sid","KEGGID")
                           }
            }else{
                
                resmat<-as.data.frame(cbind(names(kegg_ids),kegg_ids))

                colnames(resmat)<-c("pubchem_sid","KEGGID")
                resmat$pubchem_sid<-gsub(resmat$pubchem_sid,pattern="pubchem:",replacement="")
                resmat$KEGGID<-gsub(resmat$KEGGID,pattern="cpd:",replacement="")
                rownames(resmat)<-NULL
                resmat<-unique(resmat)
                
                res_final<-cbind(d1[[1]]$Information[[j]]$CID,resmat)
                colnames(res_final)<-c("pubchem_cid","pubchem_sid","KEGGID")
                
                
                
            }
            
            res_final<-as.data.frame(res_final)
            return(res_final)
		}
	})


    res_final<-ldply(sid_list,rbind)
    colnames(res_final)<-c("pubchem_cid","pubchem_sid","KEGGID")
    res_final<-as.data.frame(res_final)
	return(res_final)

	}

	}

}



#input: 
#1. KEGG species code (e.g. "hsa" for homo sapiens or "rno" for rats)
#2. KEGG database name (e.g. "pathway")
#output:
#Mapping between KEGG pathway id, gene or compound ids, and names
get_kegg_gene_compound_pathway_table<-function(kegg_species_code="hsa",kegg_database="pathway"){

	path_list<-keggList(kegg_species_code,database=kegg_database)
	path_ids<-names(path_list)
	path_ids<-gsub(path_ids,pattern="path:",replacement="")


	res<-lapply(1:length(path_ids),function(j,kegg_species_code){

		Sys.sleep(0.1)
		
        	k1<-keggGet(dbentries=path_ids[j])
        
		kegg_comp_mat<-{}
		kegg_gene_mat<-{}
        	kegg_comp_list<-k1[[1]]$COMPOUND
		kegg_gene_list<-k1[[1]]$GENE
		
		if(length(kegg_comp_list)>0){
		kegg_comp_mat<-cbind("compound",names(kegg_comp_list),kegg_comp_list)
		kegg_comp_mat<-as.data.frame(kegg_comp_mat)

		
			colnames(kegg_comp_mat)<-c("id_type","id","name")
		}

		if(length(kegg_gene_list)>0){
		kegg_gene_ids<-kegg_gene_list[seq(1,length(kegg_gene_list),2)]
		kegg_gene_names<-kegg_gene_list[seq(2,length(kegg_gene_list),2)]
		
		#get NCBI Gene ID from KEGG Gene ID 							
		ncbi_geneid<-keggConv(target="ncbi-geneid",paste(kegg_species_code,":",kegg_gene_ids,sep="")) 

		
		
		#get NCBI Protein ID from KEGG Gene ID; 
		ncbi_proteinid=keggConv(target="ncbi-proteinid",paste(kegg_species_code,":",kegg_gene_ids,sep=""))		
		

		kegg_gene_ids<-paste(kegg_gene_ids,";",ncbi_geneid,";",ncbi_proteinid,sep="")
		
		kegg_gene_mat<-cbind("gene",kegg_gene_ids,kegg_gene_names)
		kegg_gene_mat<-as.data.frame(kegg_gene_mat)

		colnames(kegg_gene_mat)<-c("id_type","id","name")

		}

		kegg_resmat<-rbind(kegg_comp_mat,kegg_gene_mat)

		kegg_resmat<-as.data.frame(kegg_resmat)

		if(length(kegg_resmat)>0){
			kegg_resmat<-cbind(path_ids[j],kegg_resmat)

			colnames(kegg_resmat)<-c("kegg_pathway_id","id_type","id","name")

			return(kegg_resmat)
		}

	},kegg_species_code=kegg_species_code)
	
	res_all<-ldply(res,rbind)
	return(res_all)
}


#Downloads and merges pathway memberships using NCBI Gene IDs, ChEBI IDs, and 
get_reactome_atlas<-function(species_code="Homo sapiens"){

	gene_to_pathways<-read.delim("https://reactome.org/download/current/NCBI2Reactome_PE_All_Levels.txt",header=FALSE)
	chebi_to_pathways<-read.delim("https://reactome.org/download/current/ChEBI2Reactome_PE_All_Levels.txt",header=FALSE)

	uniprot_to_pathways<-read.delim("https://reactome.org/download/current/UniProt2Reactome_PE_All_Levels.txt",header=FALSE)

	reactome_chebi_hsa<-chebi_to_pathways[which(chebi_to_pathways[,8]==species_code),]
	reactome_gene_hsa<-gene_to_pathways[which(gene_to_pathways[,8]==species_code),]
	reactome_protein_hsa<-uniprot_to_pathways[which(uniprot_to_pathways[,8]==species_code),]

	reactome_chebi_hsa<-cbind(reactome_chebi_hsa[,1:7],"compound")
	reactome_gene_hsa<-cbind(reactome_gene_hsa[,1:7],"gene")
	reactome_protein_hsa<-cbind(reactome_protein_hsa[,1:7],"protein")

	
	
	colnames(reactome_chebi_hsa)<-c("ID","Stable_Identifier","Name","PathwayID","PathwayURL","PathwayName","Evidence","IDtype")
	colnames(reactome_gene_hsa)<-c("ID","Stable_Identifier","Name","PathwayID","PathwayURL","PathwayName","Evidence","IDtype")

	colnames(reactome_protein_hsa)<-c("ID","Stable_Identifier","Name","PathwayID","PathwayURL","PathwayName","Evidence","IDtype")


	reactome_atlas<-rbind(reactome_chebi_hsa,reactome_gene_hsa,reactome_protein_hsa)

	return(reactome_atlas)


}

#input: PubChem CIDs
#output: data frame with PubChem CIDs mapped to CHEBI IDs
pubchem_to_chebi<-function(cids){

	
	r1<-cts_convert(query=cids,from="PubChem CID",to="InChiKey")

	chebiid<-get_chebiid(query=r1,from="INCHI/INCHI KEY")
	
	chebiid<-na.omit(chebiid)
	
	return(chebiid)
}


########END FUNCTIONS########



#load MoTrPAC data dictionary
g1<-MotrpacBicQC::get_and_validate_mdd()

cids<-g1$pubchem_cid #c(g1$pubchem_cid[1:10],"6137","5283137","5281912","5280724","5281")

#map PubChem IDs to KEGG IDs
kegg_mapping_res<-pubchem_to_kegg(cids) #(g1$pubchem_cid)

#add KEGG IDs to the original table
g2=merge(g1,kegg_mapping_res,by.x="pubchem_cid",by.y="pubchem_cid",all=TRUE)

#map PubChem IDs to ChEBI IDs
chebi_mapping_res<-pubchem_to_chebi(cids) #(g1$pubchem_cid)

#add ChEBI IDs to the original table
g3=merge(g2,chebi_mapping_res[,c("chebiid","query")],by.x="inchi_key",by.y="query",all=TRUE)

g3$chebiid<-gsub(g3$chebiid,pattern="CHEBI:",replacement="")

g3<-g3[order(as.numeric(as.character(g3$chebiid))),]

g4<-g3[,c("refmet_name","super_class","main_class","sub_class","formula","exactmass","inchi_key","pubchem_cid","standard","pubchem_sid","KEGGID","chebiid")]

g4u<-g4[order(g4$pubchem_cid,g4$KEGGID),]

g4u<-g4u[-which(duplicated(paste(g4u$pubchem_cid,"_",g4u$inchi_key,sep=""))==TRUE),]

g4u2<-g4[-which(g4$refmet_name%in%g4u$refmet_name),]
g5<-rbind(g4u,g4u2)

check_nas<-which(is.na(g5$refmet_name)==TRUE)
if(length(check_nas)>0){
g5<-g5[-which(is.na(g5$refmet_name)==TRUE),]
}

write.csv(g5,file="updated_metabolite_data_dictionary_08062020.csv")


#get KEGG molecular atlas for humans and rats
kegg_atlas_hsa<-get_kegg_gene_compound_pathway_table(kegg_species_code="hsa",kegg_database="pathway")
kegg_atlas_rno<-get_kegg_gene_compound_pathway_table(kegg_species_code="rno",kegg_database="pathway")


#get Reactome molecular atlas for humans and rats
reactome_hsa<-get_reactome_atlas("Homo sapiens")
reactome_rno<-get_reactome_atlas("Rattus norvegicus")
