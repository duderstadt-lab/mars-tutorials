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

//Script B: beta/gamma correction with two populations after the archive merge


//**********************************


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
//Interpretation of the I vs T traces and E and S calculation commands
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









