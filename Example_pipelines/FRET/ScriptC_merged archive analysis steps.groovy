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

//Table of contents
//1. Add the region 'active' to each metadata archive in the merged archive
//2. Filter the molecules based on the outcome of the KCP analysis and tag accordingly
//3. Calculate the apparent intensity values from the intensity vs. T traces
//4. Calculate the uncorrected apparent E and S values (iEapp and iSapp)
//5. Apply an archive-wise trace-wise background correction
//6. Calculate the alpha and delta correction factors and apply the correction on E and S
//7. Calculate the beta and gamma correction factors and apply the correction on E and 

//******************************************************************************

//1. Add the region 'active' to each metadata archive in the merged archive

//import dependencies
#@ MoleculeArchive archive

import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.util.*

//add the region
archive.metadata().forEach{metadata ->
	metadata.putRegion(new MarsRegion("active", "T", 2, 500, "#416EF4", 0.2))
}

//2. Filter the molecules based on the outcome of the KCP analysis and tag accordingly

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


//tag the molecules according to their bleaching parameter values as calculated previously
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

//3. Calculate the apparent intensity values from the intensity vs. T traces

//import dependencies
#@ MoleculeArchive archive

import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import de.mpg.biochem.mars.util.*
import org.scijava.table.*

//determine the intensity prior to bleaching as identified in the KCP analysis
//add this value as parameter to each molecule record
archive.getMoleculeUIDs().stream().forEach({UID ->
  Molecule molecule = archive.get(UID)
  if (molecule.hasSegmentsTable("T","0","active")){
  MarsTable table = molecule.getSegmentsTable("T","0","active")
  double ymax = table.max("A")
  molecule.setParameter("Iaemaex",ymax)}
  if (molecule.hasSegmentsTable("T","1 Green","active")){
  MarsTable table2 = molecule.getSegmentsTable("T","1 Green","active")
  double ymax2 = table2.max("A")
  molecule.setParameter("Idemdex",ymax2)}
  if (molecule.hasSegmentsTable("T","1 Red","active")){
  MarsTable table3 = molecule.getSegmentsTable("T","1 Red","active")
  double ymax3 = table3.max("A")
  molecule.setParameter("Iaemdex",ymax3)}
  }) 

//4. Calculate the uncorrected apparent E and S values (iEapp and iSapp)

//import dependencies
#@ MoleculeArchive archive

import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import de.mpg.biochem.mars.util.*
import org.scijava.table.*

//retrieve the intensity parameters and calculate the iEapp and iSapp value for each molecule
//store the iEapp and iSapp values as a parameter in the respective molecule record
archive.getMoleculeUIDs().stream().forEach({UID ->
  Molecule molecule = archive.get(UID)
  if (molecule.hasTag("Active_single")){
  double Iaemdex = molecule.getParameter("Iaemdex")
  double Idemdex = molecule.getParameter("Idemdex")
  double Iaemaex = molecule.getParameter("Iaemaex")
  double S = (Iaemdex + Idemdex) / (Iaemdex + Idemdex + Iaemaex)
  double E = Iaemdex / (Iaemdex + Idemdex)
  molecule.setParameter("iSapp",S)
  molecule.setParameter("iEapp",E)}
  })


//5. Apply an archive-wise trace-wise background correction

//import dependencies
#@ MoleculeArchive archive

import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import de.mpg.biochem.mars.util.*
import org.scijava.table.*
import groovy.lang.*

//create a list containing the names of all videos present in the Archive
//make this a unique list
def names_list = []

archive.metadata().forEach{metadata ->
	names_list.add(metadata.getImage(0).getName())
}
names_list.unique()

//iterate through all videos
//for each video, create a list with the UIDs of the three corresponding archives [FRET, AO, DO]
//calculate the background values from the FRET archive and correct all archives with this value
for (def name : names_list) {

	//create the [FRET, AO, DO] UID list for the video
	String[] metadataUIDs = new String[3]
	for (def mUID1 : archive.getMetadataUIDs()) {
		def metadata = archive.getMetadata(mUID1)
		if (metadata.getImage(0).getName().equals(name)) {
				if (metadata.hasTag("FRET"))
					metadataUIDs[0] = mUID1
				else if (metadata.hasTag("AO"))
					metadataUIDs[1] = mUID1
				else if (metadata.hasTag("DO"))
					metadataUIDs[2] = mUID1
		}
	}

	//calculate the correction factors from the FRET archive
	if (metadataUIDs[0] != null) {
		def metadata = archive.getMetadata(metadataUIDs[0]) //select the FRET archive
 		def background_aemaex = []
		def background_aemdex = []
		def background_demdex = []
		archive.molecules()
		.filter{molecule -> molecule.getMetadataUID().equals(metadataUIDs[0])}
		.filter{molecule -> molecule.hasTag("Active_single")}	//only use molecules tagged 'Active_single'
		.forEach{molecule ->		//loop through all molecules
			if (molecule.hasSegmentsTable("T","0","active")) {
	  			MarsTable table = molecule.getSegmentsTable("T","0","active")
	  			double Ibackground = table.min("A")
	  			background_aemaex.add(Ibackground)
			}
	  		if (molecule.hasSegmentsTable("T","1 Green","active")) {
	  			MarsTable table2 = molecule.getSegmentsTable("T","1 Green","active")
	  			double Ibackground2 = table2.min("A")
	  			background_demdex.add(Ibackground2)
			}
	  		if (molecule.hasSegmentsTable("T","1 Red","active")) {
	  			MarsTable table3 = molecule.getSegmentsTable("T","1 Red","active")
	  			double Ibackground3 = table3.min("A")
	  			background_aemdex.add(Ibackground3)
	  	 	}
  		}
 
  	//calculate the background correction factors and store in the metadata
  	double mean_aemaex = background_aemaex.sum() / background_aemaex.size()	//calculate the mean background values
	double mean_aemdex = background_aemdex.sum() / background_aemdex.size()
	double mean_demdex = background_demdex.sum() / background_demdex.size()	
	metadata.setParameter("background_aemaex", mean_aemaex) 
	metadata.setParameter("background_aemdex", mean_aemdex)
	metadata.setParameter("background_demdex", mean_demdex)

	//add the background correction factors to the AO and DO archive as well
	if (metadataUIDs[0] != null && metadataUIDs[1] != null) {
		def metadataAO = archive.getMetadata(metadataUIDs[1]) //select the AO archive
		metadataAO.setParameter("background_aemaex", mean_aemaex) 
		metadataAO.setParameter("background_aemdex", mean_aemdex)
		metadataAO.setParameter("background_demdex", mean_demdex)
	}

	if (metadataUIDs[0] != null && metadataUIDs[2] != null) {
		def metadataDO = archive.getMetadata(metadataUIDs[2]) //select the DO archive
		metadataDO.setParameter("background_aemaex", mean_aemaex) 
		metadataDO.setParameter("background_aemdex", mean_aemdex)
		metadataDO.setParameter("background_demdex", mean_demdex)
	}	
	}
}

//correct all I values with the background correction factor stored in the metadata record
for (def mUID : archive.getMetadataUIDs()) { //iterate through metadata records
	def metadata = archive.getMetadata(mUID)
	double background_Iaemaex = metadata.getParameter("background_aemaex")
	double background_Iaemdex = metadata.getParameter("background_aemdex")
	double background_Idemdex = metadata.getParameter("background_demdex")
	println(background_Iaemaex)

	archive.molecules()
	.filter{molecule -> molecule.getMetadataUID().equals(mUID)}
	.forEach{molecule -> //iterate through molecule records
		double Iaemdex = molecule.getParameter("Iaemdex")
  		double Idemdex = molecule.getParameter("Idemdex")
 		double Iaemaex = molecule.getParameter("Iaemaex")
 		double iiIaemdex = Iaemdex - background_Iaemdex
  		double iiIdemdex = Idemdex - background_Idemdex
  		double iiIaemaex = Iaemaex - background_Iaemaex
  		molecule.setParameter("iiIaemdex",iiIaemdex)
  		molecule.setParameter("iiIdemdex",iiIdemdex)
  		molecule.setParameter("iiIaemaex",iiIaemaex)
	}
	
}

//6. Calculate the alpha and delta correction factors and apply the correction on E and S

//import dependencies
#@ MoleculeArchive archive

import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import de.mpg.biochem.mars.util.*
import org.scijava.table.*

def E_list = []
def S_list = []

//calculate the corrected S and E values for the AO and DO populations respectively
archive.getMoleculeUIDs().stream().forEach({UID ->
  def molecule = archive.get(UID)
  if (molecule.hasTag("AO_active")){
  	double Iaemdex = molecule.getParameter("iiIaemdex")
    double Idemdex = molecule.getParameter("iiIdemdex")
    double Iaemaex = molecule.getParameter("iiIaemaex")
  	double S = (Iaemdex + Idemdex) / (Iaemdex + Idemdex + Iaemaex)
  	double E = (Iaemdex) / (Iaemdex + Idemdex)
    molecule.setParameter("iiSappAO",S)
    molecule.setParameter("iiEappAO",E)
    S_list.add(S)
  	}
  	else if (molecule.hasTag("DO_active")){
		double Iaemdex = molecule.getParameter("iiIaemdex")
        double Idemdex = molecule.getParameter("iiIdemdex")
        double Iaemaex = 0
  		double E = Iaemdex / (Iaemdex + Idemdex)
  		double S = (Iaemdex + Idemdex) / (Iaemdex + Idemdex + Iaemaex)
  		molecule.setParameter("iiEappDO",E)
  		molecule.setParameter("iiSappDO",S)
  		E_list.add(E)
  	}
  })

//calculate the mean S and E values from the AO and DO populations and assign a parameter in the metadata
double E_mean = E_list.stream().filter{item -> item != Double.NaN}.mapToDouble{item -> Double.valueOf(item)}.sum()/E_list.size()
double S_mean =  S_list.stream().filter{item -> item != Double.NaN}.mapToDouble{item -> Double.valueOf(item)}.sum()/S_list.size()

archive.metadata().forEach{metadata -> 
	metadata.setParameter("iiEappDO_mean", E_mean)
	metadata.setParameter("iiSappAO_mean", S_mean)
}


//calculate alpha and delta and assign a parameter in the metadata
double alpha = E_mean/(1- E_mean)
double delta = S_mean/(1-S_mean)

archive.metadata().forEach{metadata -> 
	metadata.setParameter("alpha", alpha)
	metadata.setParameter("delta", delta)
}

//calculate FAD for each FRET molecule and the corrected E and S values
archive.getMoleculeUIDs().stream().forEach({UID ->
  def molecule = archive.get(UID)
  if (molecule.hasTag("Active_single")){
  	double Iaemdex = molecule.getParameter("iiIaemdex")
    double Idemdex = molecule.getParameter("iiIdemdex")
    double Iaemaex = molecule.getParameter("iiIaemaex")
    double FAD = Iaemdex - alpha * Idemdex - delta * Iaemaex
    double iiiE = FAD / (FAD + Idemdex)
    double iiiS = (FAD + Idemdex) / (FAD + Idemdex + Iaemaex)
    molecule.setParameter("F_AD",FAD)
    molecule.setParameter("iiiEapp",iiiE)
    molecule.setParameter("iiiSapp",iiiS)
  }
})

//7. Calculate the beta and gamma correction factors and apply the correction on E and S

//import dependencies
#@ MoleculeArchive archive

import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import de.mpg.biochem.mars.util.*
import org.scijava.table.*
import org.apache.commons.math3.stat.regression.SimpleRegression

def regression = new SimpleRegression(true)
def E_list = []
def S_list = []

//create a list with all data points (iiiEapp and iiiSapp) to be used for the fit
archive.getMoleculeUIDs().stream().forEach({UID ->
  def molecule = archive.get(UID)
  if (molecule.hasTag("Active_single")){
  	double E = molecule.getParameter("iiiEapp")
    double S = 1 / molecule.getParameter("iiiSapp")
    S_list.add(S)
    E_list.add(E)}})

//add the datapoints to be included in the regression analysis and store found parameters
for (i in [E_list,S_list].transpose()){
	regression.addData(i[0],i[1])
}

double a = regression.getSlope()
double b = regression.getIntercept()
double beta = a + b - 1
double gamma = (b - 1) / (a + b - 1)

archive.metadata().forEach{metadata -> 
	metadata.setParameter("beta", beta)
	metadata.setParameter("gamma", gamma)
}

//calculate F_AA and F_DD values for each molecule and calculate the E and S values
archive.getMoleculeUIDs().stream().forEach({UID ->
  def molecule = archive.get(UID)
  if (molecule.hasTag("Active_single")){
  	double FAD = molecule.getParameter("F_AD")
  	double Idemdex = molecule.getParameter("iiIdemdex")
  	double Iaemaex = molecule.getParameter("iiIaemaex")
  	double FDD = gamma * Idemdex
  	double FAA = Iaemaex / beta
    molecule.setParameter("F_DD",FDD)
    molecule.setParameter("F_AA",FAA)
    double E = FAD / (FDD + FAD)
    double S = (FDD + FAD) / (FDD + FAD + FAA)
    molecule.setParameter("E",E)
    molecule.setParameter("S",S)
  }
    
})