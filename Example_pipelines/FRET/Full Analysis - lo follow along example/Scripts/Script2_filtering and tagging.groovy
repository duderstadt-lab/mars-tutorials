/*******************************************************************************
 * Copyright (C) 2021, Duderstadt Lab
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 ******************************************************************************/
//written by: Nadia M. Huisjes, MSc. (Duderstadt lab)

//This script accompanies the 'FRET dataset analysis using Mars' example pipeline as described on the mars docs.
//https://duderstadt-lab.github.io/mars-docs/examples/FRET

// Filter the molecules based on the outcome of the KCP analysis and tag accordingly

//import dependencies
#@ MoleculeArchive archive

import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import de.mpg.biochem.mars.util.*
import org.scijava.table.*

//asses whether the red dye bleached in the trace and assign Boolean parameters
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

//asses whether the green dye bleached in the trace and assign Boolean parameters
//discriminate between two options: A bleach before D, results in 3 D states
//									A bleach after D, results in 2 D states


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
  		else{molecule.setParameter("Green_step", false)}
  	if (Seg_diff>0 && steps==2){
		molecule.setParameter("Green_dual_state", true)}
		else{molecule.setParameter("Green_dual_state", false)}
	if (steps==3){
		double ymiddle = table.getColumnAsDoubles("A")[1]
		if (ymin<ymax && ymax<ymiddle){
			molecule.setParameter("Green_tri_state", true)
		}	else{molecule.setParameter("Green_tri_state", false)}
		}
		else{molecule.setParameter("Green_tri_state", false)}
  }
  	else{
  		molecule.setParameter("Green_dual_state", false)
  		molecule.setParameter("Green_tri_state", false)
  		molecule.setParameter("Green_step",false)}
})



//asses whether the bleaching was observed in the FRET trace and assign Boolean parameters
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

//calculate the difference between the green and red signal after photobleaching
//these signals should be relatively similar to one another, else the peak should be excluded from the analysis.
archive.getMoleculeUIDs().stream().forEach({UID ->
  Molecule molecule = archive.get(UID)
  if (molecule.hasSegmentsTable("T","0","active") && molecule.hasSegmentsTable("T","1 Green","active")){
  	MarsTable table = molecule.getSegmentsTable("T","0","active")
  	MarsTable table2 = molecule.getSegmentsTable("T","1 Green","active")
  	double Ibckaemaex = table.min("A")
  	double Ibckdemdex = table2.min("A")
  	double factor = Math.abs(Ibckaemaex / Ibckdemdex)
  	double factor2 = Math.abs(Ibckdemdex / Ibckaemaex)
  	if (factor > 2 || factor2 > 2){
  		molecule.addTag("Outlier_signal")}
  	}
  }
)


//tag the molecules according to their bleaching parameter values as calculated previously
archive.getMoleculeUIDs().stream().forEach({UID ->
  def molecule = archive.get(UID)
  if (!molecule.hasTag("Outlier_signal")){
  	if (archive.metadataHasTag(molecule.getMetadataUID(),"FRET")){
  	 if (molecule.getBooleanParameter("Red_singlebleach") && molecule.getBooleanParameter("FRET_bleach") && molecule.getBooleanParameter("Green_dual_state")){
  		 molecule.addTag("Active_single")}
  	 if (molecule.getBooleanParameter("Red_singlebleach") && molecule.getBooleanParameter("FRET_bleach") && molecule.getBooleanParameter("Green_tri_state")){
  		 molecule.addTag("Active_single")}
   }
   if (archive.metadataHasTag(molecule.getMetadataUID(),"AO")){
  	 if (molecule.getBooleanParameter("Red_singlebleach")&& !molecule.getBooleanParameter("Green_step")&& !molecule.getBooleanParameter("FRET_bleach")){
  		 molecule.addTag("AO_active")}}
   if (archive.metadataHasTag(molecule.getMetadataUID(),"DO")){
  	 if (molecule.getBooleanParameter("Green_dual_state")&& !molecule.getBooleanParameter("Red_step") && !molecule.getBooleanParameter("FRET_bleach")){
  		 molecule.addTag("DO_active")}
  		}
  }
})
