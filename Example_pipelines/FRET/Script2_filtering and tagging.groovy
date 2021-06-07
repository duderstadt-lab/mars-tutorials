#@ MoleculeArchive archive

import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import de.mpg.biochem.mars.util.*
import org.scijava.table.*

//Asses whether the red dye bleached in the trace and assign Boolean parameters
archive.getMoleculeUIDs().stream().forEach({UID ->
  Molecule molecule = archive.get(UID)
  if (molecule.hasSegmentsTable("T","0","active")){
	MarsTable table = molecule.getSegmentsTable("T","0","active")
  	double ymax = table.getColumnAsDoubles("A")[0]
	double ymin = table.getColumnAsDoubles("A")[-1] 	
  	double Seg_diff = ymax - ymin
  	double steps = table.getColumnAsDoubles("A").size()
  	if (Seg_diff != 0){
  	molecule.setParameter("Red_step",true)}
  		else{
  			molecule.setParameter("Red_step",false)}
  	if (Seg_diff>0 && steps==2){
  		molecule.setParameter("Red_singlebleach", true)}
  		else{molecule.setParameter("Red_singlebleach", false)}
  }
  	else{
  		molecule.setParameter("Red_singlebleach", false)
  		molecule.setParameter("Red_step",false)}
})

//Asses whether the green dye bleached in the trace and assign Boolean parameters
archive.getMoleculeUIDs().stream().forEach({UID ->
  Molecule molecule = archive.get(UID)
  if (molecule.hasSegmentsTable("T","1 Green","active")){
	MarsTable table = molecule.getSegmentsTable("T","1 Green","active")
  	double ymax = table.getColumnAsDoubles("A")[0]
	double ymin = table.getColumnAsDoubles("A")[-1]
  	double Seg_diff = ymax - ymin
  	double steps = table.getColumnAsDoubles("A").size()
  	if (Seg_diff !=0){
  	molecule.setParameter("Green_step",true)}
  		else{
  			molecule.setParameter("Green_step",false)}
  	if (Seg_diff>0 && steps==2){
  		molecule.setParameter("Green_singlebleach", true)}
  		else{molecule.setParameter("Green_singlebleach", false)}
  }
  	else{
  		molecule.setParameter("Green_singlebleach", false)
  		molecule.setParameter("Green_step",false)}
})

//Asses whether the bleaching was observed in the FRET trace and assign Boolean parameters

archive.getMoleculeUIDs().stream().forEach({UID ->
  Molecule molecule = archive.get(UID)
  if (molecule.hasSegmentsTable("T","1 Red","active")){
  	MarsTable table = molecule.getSegmentsTable("T","1 Red","active")
  	double ymin = table.min("A")
  	double ymax = table.max("A")
  	double diff = ymax - ymin
  	if (diff==0){
  		molecule.setParameter("FRET_bleach",false)
  	} else{
  		molecule.setParameter("FRET_bleach",true)}}
})


//Tag the molecules according to their bleaching parameter values as calculated previously
archive.getMoleculeUIDs().stream().forEach({UID ->
  def molecule = archive.get(UID)
  if (archive.metadataHasTag(molecule.getMetadataUID(),"FRET")){
  	if (molecule.getBooleanParameter("Red_step")&& molecule.getBooleanParameter("Green_step")&& molecule.getBooleanParameter("FRET_bleach")){
  	molecule.addTag("Active")
    } else {
  	molecule.addTag("background")}
  if (molecule.hasTag("Active")&&molecule.getBooleanParameter("Red_singlebleach")&&molecule.getBooleanParameter("Green_singlebleach")&& molecule.getBooleanParameter("FRET_bleach")){
  	molecule.addTag("Active_single")}}
  if (archive.metadataHasTag(molecule.getMetadataUID(),"AO")){
  	if (molecule.getBooleanParameter("Red_singlebleach")&& !molecule.getBooleanParameter("Green_step")&& !molecule.getBooleanParameter("FRET_bleach")){
  		molecule.addTag("AO_active")}}
  if (archive.metadataHasTag(molecule.getMetadataUID(),"DO")){
  	if (molecule.getBooleanParameter("Green_singlebleach")&& !molecule.getBooleanParameter("Red_step") && !molecule.getBooleanParameter("FRET_bleach")){
  		molecule.addTag("DO_active")}}
})

