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



//**********************************
//------User Settings-------//

//Set archive type
molecule_kind = "mid" 	//"lo" or "mid" ; determines which metadata tag is set
videos = "0-6"			//"0-4" or "5-6", "0-6" ; determines which KCP settings are used
filter_S_tres = 20
filter_E_tres = 20
filter_bool = true		// filter values before mean/median calculation for alpha/delta
median_bool = true		// or false (then the mean will be used for the alpha/delta calculation)


//**********************************
//No user changes needed beyond this point


//Import packages
import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import de.mpg.biochem.mars.util.*
import org.scijava.table.*
import de.mpg.biochem.mars.kcp.commands.*
import groovy.lang.*
import org.apache.commons.math3.stat.regression.SimpleRegression


//Interaction parameters
#@ MoleculeArchive archive
#@ ImageJ ij


//------Method definitions-------//
//Add metadata tag command
def MetadataTag(String name){
	archive.metadata().forEach{metadata ->
		metadata.addTag(name)
	}
}

//KCP command
def RunKCP(String KCP_Y, BigDecimal KCP_conf, Integer KCP_sigma, String KCP_region){
	final KCPCommand kcpCalc = new KCPCommand()
	kcpCalc.setContext(ij.getContext())
	kcpCalc.setArchive(archive)
	kcpCalc.setXColumn("T")
	kcpCalc.setYColumn(KCP_Y)
	kcpCalc.setConfidenceLevel(KCP_conf)
	kcpCalc.setGlobalSigma(KCP_sigma)
	kcpCalc.setRegionSource("Metadata")
	kcpCalc.setCalculateBackgroundSigma(false)
	kcpCalc.setBackgroundRegion("None")
	kcpCalc.setAnalyzeRegion(true)
	kcpCalc.setRegion(KCP_region)
	kcpCalc.setFitSteps(true)
	kcpCalc.setIncludeTags("All")
	kcpCalc.setTags("None")
	kcpCalc.setThreads(8)
	kcpCalc.run()
}

//Segment Table based filtering condition calculation
def Seg_diff(def molecule, String Ycol){
	MarsTable table = molecule.getSegmentsTable("T",Ycol,"active")
  	double ymax = table.getColumnAsDoubles("A")[0]
	double ymin = table.getColumnAsDoubles("A")[-1]
	double diff = ymax - ymin 	
  	return diff
}

def Seg_steps(def molecule, String Ycol){
	MarsTable table = molecule.getSegmentsTable("T",Ycol,"active")
	double seg_steps = table.getColumnAsDoubles("A").size()
	return seg_steps
}

def Seg_middle_diff(def molecule, String Ycol){
	MarsTable table = molecule.getSegmentsTable("T",Ycol,"active")
	double seg_middle = table.getColumnAsDoubles("A")[1]
	double ymax = table.getColumnAsDoubles("A")[0]
	double seg_middle_diff = seg_middle - ymax 
	return seg_middle_diff
}

//Interpretation of the I vs T traces and E and S calculation commands

def Find_T_endfret(def molecule, String Ycol){
	MarsTable table = molecule.getSegmentsTable("T",Ycol,"active")
	double T = table.getValue("X1",1)
	return T
}

def Calc_E_S(String E_name, String S_name, String Iaemdex, String Idemdex, String Iaemaex){
	archive.getMoleculeUIDs().stream().forEach({UID ->
  	Molecule molecule = archive.get(UID)
  	if (molecule.hasTag("Active_single") || molecule.hasTag("AO_active") || molecule.hasTag("DO_active")){
  		double Tendfret = molecule.getParameter("Tendfret")
  		MarsTable table = molecule.getTable()
  		double len = table.getRowCount()
  		for (i=0; i<Tendfret; i++){
    		double Iaemaex_val = table.getValue(Iaemaex,i)
    		double Iaemdex_val = table.getValue(Iaemdex,i)
    		double Idemdex_val = table.getValue(Idemdex,i)
  			double S = (Iaemdex_val + Idemdex_val) / (Iaemdex_val + Idemdex_val + Iaemaex_val)
    		double E = Iaemdex_val / (Iaemdex_val + Idemdex_val)
    		table.setValue(E_name,i,E)
    		table.setValue(S_name,i,S)
  		}
  		for (i=Tendfret; i<len; i++){
  		table.setValue(E_name,(int)i,"NaN")
    	table.setValue(S_name,(int)i,"NaN")
  		}
  	}
  	})
}

//Methods for background correction of the traces
def VideoList (){
	//Returns a unique list with video file names
	def names_list = []
	archive.metadata().forEach{metadata ->
	names_list.add(metadata.getImage(0).getName())
	}
	return names_list.unique()
}

def BackgroundIntensity(def molecule, String Ycol){
	MarsTable table = molecule.getSegmentsTable("T",Ycol,"active")
    double Ibackground = table.min("A")
    return Ibackground
}

def Mean(List l1){
	double sum = l1.sum()
	double len = l1.size()
	double mean = sum / len
	return mean
}

def Median(List l1){
	double len = l1.size()
	double modulo = len % 2
	double ind = Math.floor(len/2)
	if (modulo == 0){
		double var1 = l1[ind]
		double var2 = l1[ind - 1]
		double median = (var1 + var2) / 2
		return median
	} else{
		double median = l1[ind]
		return median
	}
}

def Background_Video_based(){
	names_list = VideoList()
	for (def name : names_list) {
	String[] metadataUIDs = new String[3]
	for (def mUID1 : archive.getMetadataUIDs()) {
		def metadata = archive.getMetadata(mUID1)
				if (metadata.getImage(0).getName().equals(name)) {
				if (metadata.hasTag("FRET")){
					metadataUIDs[0] = mUID1}
				else if (metadata.hasTag("AO")){
					metadataUIDs[1] = mUID1}
				else if (metadata.hasTag("DO")){
					metadataUIDs[2] = mUID1}
		}
	}
	if (metadataUIDs[0] != null){
		def metadata = archive.getMetadata(metadataUIDs[0])
		def background_aemaex = []
		def background_aemdex = []
		def background_demdex = []
		archive.molecules()
		.filter{molecule -> molecule.getMetadataUID().equals(metadataUIDs[0])}
		.filter{molecule -> molecule.hasTag("Active_single")}
		.forEach{molecule ->
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
	
	double mean_aemaex = Mean(background_aemaex)
	double mean_aemdex = Mean(background_aemdex)
	double mean_demdex = Mean(background_demdex)
	metadata.setParameter("background_aemaex", mean_aemaex) 
	metadata.setParameter("background_aemdex", mean_aemdex)
	metadata.setParameter("background_demdex", mean_demdex)
	
	//add the mean values to the AO archive
	if (metadataUIDs[0] != null && metadataUIDs[1] != null) {
	def metadataAO = archive.getMetadata(metadataUIDs[1]) 
	metadataAO.setParameter("background_aemaex", mean_aemaex) 
	metadataAO.setParameter("background_aemdex", mean_aemdex)
	metadataAO.setParameter("background_demdex", mean_demdex)
	}
	
	//add the mean values to the DO archive
	if (metadataUIDs[0] != null && metadataUIDs[2] != null) {
	def metadataDO = archive.getMetadata(metadataUIDs[2]) 
	metadataDO.setParameter("background_aemaex", mean_aemaex) 
	metadataDO.setParameter("background_aemdex", mean_aemdex)
	metadataDO.setParameter("background_demdex", mean_demdex)
	}
	}
	}
}


def I_background_corr_metadata_based(){
	for (def mUID : archive.getMetadataUIDs()) { 
	def metadata = archive.getMetadata(mUID)
	double background_Iaemaex = metadata.getParameter("background_aemaex")
	double background_Iaemdex = metadata.getParameter("background_aemdex")
	double background_Idemdex = metadata.getParameter("background_demdex")

	archive.molecules()
	.filter{molecule -> molecule.getMetadataUID().equals(mUID)}
	.forEach{molecule -> 
    if (molecule.hasTag("Active_single") || molecule.hasTag("AO_active") || molecule.hasTag("DO_active")){
      double Tendfret = molecule.getParameter("Tendfret")
      MarsTable table = molecule.getTable()
      double len = table.getRowCount()
      for (i=0; i<Tendfret; i++){
        double Iaemaex = table.getValue("Iaemaex",i)
        double Iaemdex = table.getValue("Iaemdex",i)
        double Idemdex = table.getValue("Idemdex",i)
        double iiIaemaex = Iaemaex - background_Iaemaex
        double iiIdemdex = Idemdex - background_Idemdex
        double iiIaemdex = Iaemdex - background_Iaemdex
        table.setValue("iiIaemaex",i,iiIaemaex)
        table.setValue("iiIdemdex",i,iiIdemdex)
        table.setValue("iiIaemdex",i,iiIaemdex)}
      for (i=Tendfret; i<len; i++){
  		table.setValue("iiIaemaex",(int)i,"NaN")
        table.setValue("iiIdemdex",(int)i,"NaN")
        table.setValue("iiIaemdex",(int)i,"NaN")}
	}
}}
}

def alpha_delta_lists(){
	//Generate the iiEapp and iiSapp list from all DO and AO molecules respectively
	def E_list = []
	def S_list = []

	archive.getMoleculeUIDs().stream().forEach({UID ->
  	Molecule molecule = archive.get(UID)
  	if (molecule.hasTag("AO_active") || molecule.hasTag("DO_active")){
  		double Tendfret = molecule.getParameter("Tendfret")
  		MarsTable table = molecule.getTable()
		for (i=0; i<Tendfret; i++){
			double E = table.getValue("iiEapp", i)
			double S = table.getValue("iiSapp", i)
			if (molecule.hasTag("AO_active")){
				S_list.add(S)
			} else if(molecule.hasTag("DO_active")){
				E_list.add(E)
			}
		}
	}
	})
	return [E_list, S_list]
}

def alpha_delta_filter(int S_tres, int E_tres){
	//Filter the iiEapp and iiSapp list by a set criterium to remove outliers
	def E_list_corr = []
	def S_list_corr = []
	E_list = alpha_delta_lists()[0]
	S_list = alpha_delta_lists()[1]

	for (i=0; i<S_list.size();i++){
		if (S_list[i]<S_tres && S_list[i]>-S_tres){
			S_list_corr.add(S_list[i])
		}
	}

	for (i=0;i<E_list.size();i++){
		if (E_list[i]<E_tres && E_list[i]>-E_tres){	
			E_list_corr.add(E_list[i])
		}
	}

	return [E_list_corr, S_list_corr]
	
}

def alpha_delta_calculation(Boolean filter, int S_tres, int E_tres, Boolean median){
	if (filter == true){
		E_list = alpha_delta_filter(S_tres,E_tres)[0]
		S_list = alpha_delta_filter(S_tres,E_tres)[1]
	} else {
		E_list = alpha_delta_lists()[0]
		S_list = alpha_delta_lists()[1]
	}

	if (median == true){
		E_stat = Median(E_list)
		S_stat = Median(S_list)
		double alpha = E_stat/(1-E_stat)
		double delta = S_stat/(1-S_stat)
		archive.metadata().forEach{metadata -> 
		metadata.setParameter("iiEappDO_median", E_stat)
		metadata.setParameter("iiSappAO_median", S_stat)
		metadata.setParameter("alpha", alpha)
		metadata.setParameter("delta", delta)
		}
		return [alpha, delta]
		
	} else{
		E_stat = Mean(E_list)
		S_stat = Mean(S_list)
		double alpha = E_stat/(1-E_stat)
		double delta = S_stat/(1-S_stat)
		
		archive.metadata().forEach{metadata -> 
		metadata.setParameter("iiEappDO_mean", E_stat)
		metadata.setParameter("iiSappAO_mean", S_stat)
		metadata.setParameter("alpha", alpha)
		metadata.setParameter("delta", delta)
		}
		return [alpha, delta]
		
	}
}


def beta_gamma_lists(){
	//Generate the iiEapp and iiSapp list from all DO and AO molecules respectively
	def E_list = []
	def S_list = []

	archive.getMoleculeUIDs().stream().forEach({UID ->
  	Molecule molecule = archive.get(UID)
  	if (molecule.hasTag("Active_single")){
  		double Tendfret = molecule.getParameter("Tendfret")
  		MarsTable table = molecule.getTable()
		for (i=0; i<Tendfret; i++){
			double E = table.getValue("iiiEapp", i)
			double S = table.getValue("iiiSapp", i)
			E_list.add(E)
			S_list.add(S)
			}
	}
	})
	return [E_list, S_list]
}


def beta_gamma_calculation(){
	def lists = beta_gamma_lists()
	def E_list = lists[0]
	def S_list = lists[1]
	//initiate the regression
	def regression = new SimpleRegression(true)
	//add data to the model
	for (i in [E_list,S_list].transpose()){
		regression.addData(i[0],1/i[1])
	}
	//save the regression outputs
	double a = regression.getSlope()
	double b = regression.getIntercept()
	double beta = a + b - 1
	double gamma = (b - 1) / (a + b - 1)

	archive.metadata().forEach{metadata -> 
		metadata.setParameter("beta", beta)
		metadata.setParameter("gamma", gamma)
		}
	return [beta, gamma]
	
}

//*******************************************
//------Code execution-------//

//1. Set region, add metadata tag and run KCP
archive.metadata().forEach{metadata ->
	metadata.putRegion(new MarsRegion("active", "T", 2, 500, "#416EF4", 0.2))
}

//comment out code not applicable to the current archive
MetadataTag(molecule_kind)

if (videos="0-4"){
	RunKCP("1 Red", 0.99, 70000, "active")
	RunKCP("1 Green", 0.99, 70000, "active")
	RunKCP("0", 0.99, 70000, "active")
} 
if (videos="5-6"){
	RunKCP("1 Red", 0.99, 30000, "active")
	RunKCP("1 Green", 0.99, 30000, "active")
	RunKCP("0", 0.99, 30000, "active")
} 
if (molecule_kind = "mid"){
	RunKCP("1 Red", 0.99, 5000, "active")
	RunKCP("1 Green", 0.99, 5000, "active")
	RunKCP("0", 0.99, 8000, "active")
} else{
	print("video input incorrect")
}


//2. Filter all traces by KCP-identified bleaching steps
//2.A Investigate the number of bleaching steps observed for all intensity trace colors
//    Assign boolean parameters based on the outcome of this analysis

archive.getMoleculeUIDs().stream().forEach({UID ->
  Molecule molecule = archive.get(UID)
  //Filter all red dye bleaches
  if (molecule.hasSegmentsTable("T","0","active")){
  	double seg_diff = Seg_diff(molecule, "0")
  	double steps = Seg_steps(molecule, "0")
  	if (seg_diff != 0){
  		molecule.setParameter("Red_step",true)}
  		else{
  			molecule.setParameter("Red_step",false)}
  	if (seg_diff>0 && steps==2){
  		molecule.setParameter("Red_singlebleach", true)}
  		else{molecule.setParameter("Red_singlebleach", false)}
  }
  	else{
  		molecule.setParameter("Red_singlebleach", false)
  		molecule.setParameter("Red_step",false)}
  //Filter all green dye bleaches
  if (molecule.hasSegmentsTable("T","1 Green","active")){
  	double seg_diff = Seg_diff(molecule, "1 Green")
  	double steps = Seg_steps(molecule, "1 Green")
  	
	if (seg_diff != 0){
		molecule.setParameter("Green_step",true)}
  		else{molecule.setParameter("Green_step", false)}
  	if (seg_diff>0 && steps==2){
		molecule.setParameter("Green_dual_state", true)}
		else{molecule.setParameter("Green_dual_state", false)}
	if (steps==3){
		double middle_diff = Seg_middle_diff(molecule, "1 Green")
		if (seg_diff>0 && middle_diff>0){
			molecule.setParameter("Green_tri_state", true)
		} else{molecule.setParameter("Green_tri_state", false)
			}
	} else {
		molecule.setParameter("Green_tri_state", false)
	}}
	else{
		molecule.setParameter("Green_dual_state", false)
  		molecule.setParameter("Green_tri_state", false)
  		molecule.setParameter("Green_step",false)}
  //Filter all observed FRET losses
  if (molecule.hasSegmentsTable("T","1 Red","active")){
  	double seg_diff = Seg_diff(molecule, "1 Red")
  	if (seg_diff == 0){
  		molecule.setParameter("FRET_bleach",false)
  	} else{
  		molecule.setParameter("FRET_bleach",true)}
  	}

})

//2.B Identify outlier artefacts by comparing the green and red signal after bleaching for each molecule
//    Add the tag "Outlier_signal" in case the ratio between the two signals is outside of the threshold.
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
})

//2.C Tag molecules according to the value of sets of bleaching parameter values 
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

//3. Extract time-wise Iapp values and deposit those in the Molecule data tables
//3.A Define the parameter Tendfret denoting the timepoint at which the FRET signal was lost (according to KCP)
archive.getMoleculeUIDs().stream().forEach({UID ->
  Molecule molecule = archive.get(UID)
  //Set T_endfret for all molecules tagged "Active_single"
  if (molecule.hasSegmentsTable("T","1 Green","active") && molecule.hasTag("Active_single")){
  	double T_endfret_one = Find_T_endfret(molecule, "1 Green")
  	double T_endfret_two = Find_T_endfret(molecule, "0")
	if (T_endfret_one>T_endfret_two){
		molecule.setParameter("Tendfret",T_endfret_two)
  	} else {
  		molecule.setParameter("Tendfret",T_endfret_one)
  }
	}
  //Set T_endfret for all molecules tagged "AO_active"
  if (molecule.hasSegmentsTable("T","1 Green","active") && molecule.hasTag("AO_active")){
  	double T_endfret = Find_T_endfret(molecule, "0")
  	molecule.setParameter("Tendfret",T_endfret)
  }
  //Set T_endfret for all molecules tagged "DO_active"
  if (molecule.hasSegmentsTable("T","1 Green","active") && molecule.hasTag("DO_active")){
  	double T_endfret = Find_T_endfret(molecule, "1 Green")
  	molecule.setParameter("Tendfret",T_endfret)
  }

})

//3.B Define the columns Iaemaex, Idemdex & Iaemdex and add the corresponding values
//    Values are only inserted until the Tendfret timepoint is reached. Afterwards NaN values are supplied.
archive.getMoleculeUIDs().stream().forEach({UID ->
  Molecule molecule = archive.get(UID)
  if (molecule.hasTag("Active_single") || molecule.hasTag("AO_active") || molecule.hasTag("DO_active")){
  double Tendfret = molecule.getParameter("Tendfret")
  MarsTable table = molecule.getTable()
  double len = table.getRowCount()
  for (i=0; i<Tendfret; i++){
	double Iaemaex = table.getValue("0",i)
	double Iaemdex = table.getValue("1 Red",i)
	double Idemdex = table.getValue("1 Green",i)
	table.setValue("Iaemaex",i,Iaemaex)
	table.setValue("Idemdex",i,Idemdex)
	table.setValue("Iaemdex",i,Iaemdex)
  }
  for (i=Tendfret; i<len; i++){
  	table.setValue("Iaemaex",(int)i,"NaN")
    table.setValue("Idemdex",(int)i,"NaN")
    table.setValue("Iaemdex",(int)i,"NaN")
  }
  }
})

//4. Calculation of iEapp and iSapp values 
Calc_E_S("iEapp", "iSapp", "Iaemdex", "Idemdex", "Iaemaex")

//5. Background correction of all intensity values
//5.A Calculate the correction factors
Background_Video_based() //one mean background intensity per video, per excitation/emission combination

//develop another method to do it molecule based instead

//5.B Correct the intensity values to obtain iiI values
I_background_corr_metadata_based()

//5.C Calculation of iiEapp and iiSapp values
Calc_E_S("iiEapp", "iiSapp", "iiIaemdex", "iiIdemdex", "iiIaemaex")

//6. Data Corrections alpha, beta, gamma, delta
//6.A Calculation of the alpha and delta factors
alpha_delta_list = alpha_delta_calculation(filter_bool, filter_S_tres, filter_E_tres, median_bool) //(filter, S_tres, E_tres, median) mean taken if median=false
double alpha = alpha_delta_list[0]
double delta = alpha_delta_list[1]


//6.B Calculation of FAD ('iiiIaemdex') values corrected for alpha and delta 
archive.getMoleculeUIDs().stream().forEach({UID ->
  Molecule molecule = archive.get(UID)
  if (molecule.hasTag("Active_single") || molecule.hasTag("AO_active") || molecule.hasTag("DO_active")){
  double Tendfret = molecule.getParameter("Tendfret")
  MarsTable table = molecule.getTable()
  double len = table.getRowCount()
  for (i=0; i<Tendfret; i++){
    double Iaemaex = table.getValue("iiIaemaex",i)
    double Iaemdex = table.getValue("iiIaemdex",i)
    double Idemdex = table.getValue("iiIdemdex",i)
  	double FAD = Iaemdex - alpha * Idemdex - delta * Iaemaex
  	table.setValue("FAD",i,FAD)
  	}
  for (i=Tendfret; i<len; i++){
  	table.setValue("FAD",(int)i,"NaN")
    }
  }
})

//6.C Calculation of iiiSapp and iiiEapp (corrected for alpha and delta)
Calc_E_S("iiiEapp", "iiiSapp", "FAD", "iiIdemdex", "iiIaemaex")

//6.D Calculation of the beta and gamma factors 
def beta_gamma = beta_gamma_calculation()
double beta = beta_gamma[0]
double gamma = beta_gamma[1]

//6.E Calculation of FAA ("iiiIaemaex") and FDD ("iiiIdemdex")
archive.getMoleculeUIDs().stream().forEach({UID ->
  Molecule molecule = archive.get(UID)
  if (molecule.hasTag("Active_single") || molecule.hasTag("AO_active") || molecule.hasTag("DO_active")){
  double Tendfret = molecule.getParameter("Tendfret")
  MarsTable table = molecule.getTable()
  double len = table.getRowCount()
  for (i=0; i<Tendfret; i++){
    double Iaemaex = table.getValue("iiIaemaex",i)
    double Idemdex = table.getValue("iiIdemdex",i)
  	double FDD = gamma * Idemdex
  	double FAA = Iaemaex / beta
    table.setValue("FDD",(int)i,FDD)
    table.setValue("FAA",(int)i,FAA)
  }
  for (i=Tendfret; i<len; i++){
  	table.setValue("FDD",(int)i,"NaN")
    table.setValue("FAA",(int)i,"NaN")
  	}
  }
})

//6.F Calculation of the fully corrected S and E values 
Calc_E_S("E", "S", "FAD", "FDD", "FAA")

////









