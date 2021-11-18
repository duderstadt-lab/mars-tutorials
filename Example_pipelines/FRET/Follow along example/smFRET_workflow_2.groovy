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


#@ MoleculeArchive archive

//Import packages
import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import de.mpg.biochem.mars.util.*
import org.scijava.table.*
import de.mpg.biochem.mars.kcp.commands.*
import groovy.lang.*
import org.apache.commons.math3.stat.regression.SimpleRegression

builder = new LogBuilder()
String log = LogBuilder.buildTitleBlock("smFRET workflow 2")
builder.addParameter("Workflow version", "0.3")
builder.addParameter("Accepted FRET", archive.molecules().filter{ m -> m.hasTag("FRET") && m.hasTag("Accepted")}.count())

//6.D Calculation of the beta and gamma factors
def beta_gamma = beta_gamma_calculation()
double beta = beta_gamma[0]
double gamma = beta_gamma[1]

//6.E Calculation of FAA ("iiiIaemaex") and FDD ("iiiIdemdex")
archive.molecules().forEach{molecule ->
    MarsTable table = molecule.getTable()
    double Tendfret = 0
  		if (molecule.hasTag("DO"))
  			Tendfret = molecule.getPosition("Donor_Bleach").getPosition()
		else if (molecule.hasTag("AO"))
			Tendfret = molecule.getPosition("Acceptor_Bleach").getPosition()
		else if (molecule.hasTag("FRET"))
			Tendfret = getTendFRET(molecule)
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

//6.F Calculation of the fully corrected S and E values
Calc_E_S("E", "S", "FAD", "FDD", "FAA")

builder.addParameter("beta", beta)
builder.addParameter("gamma", gamma)
log += builder.buildParameterList()
log += "\n" + LogBuilder.endBlock(true) + "\n"
archive.logln(log)

//Utility methods
def getTendFRET(def molecule) {
	double donorBleachT = 0
	double acceptorBleachT = 0
	if (molecule.hasPosition("Donor_Bleach"))
		donorBleachT = molecule.getPosition("Donor_Bleach").getPosition()
	if (molecule.hasPosition("Acceptor_Bleach"))
		acceptorBleachT = molecule.getPosition("Acceptor_Bleach").getPosition()

	//return (donorBleachT <= acceptorBleachT) ? donorBleachT : acceptorBleachT
	return (donorBleachT <= acceptorBleachT) ? acceptorBleachT : donorBleachT
}

def Calc_E_S(String E_name, String S_name, String Iaemdex, String Idemdex, String Iaemaex) {
	double global_E = 0
	double global_S = 0
	int global_count = 0
	archive.molecules().forEach{molecule ->
		MarsTable table = molecule.getTable()
  		double endT = 0
  		if (molecule.hasTag("DO"))
  			endT = molecule.getPosition("Donor_Bleach").getPosition()
		else if (molecule.hasTag("AO"))
			endT = molecule.getPosition("Acceptor_Bleach").getPosition()
		else if (molecule.hasTag("FRET"))
			endT = getTendFRET(molecule)

  		double mol_E = 0
  		double mol_S = 0
		int mol_obs = 0
  		for (i=0; i < endT; i++) {
    		double Iaemaex_val = table.getValue(Iaemaex,i)
    		double Iaemdex_val = table.getValue(Iaemdex,i)
    		double Idemdex_val = table.getValue(Idemdex,i)
  			double S = (Iaemdex_val + Idemdex_val) / (Iaemdex_val + Idemdex_val + Iaemaex_val)
    		double E = Iaemdex_val / (Iaemdex_val + Idemdex_val)
    		table.setValue(E_name,i,E)
    		table.setValue(S_name,i,S)
    		mol_E += E
    		mol_S += S
			mol_obs++
  		}

		if (molecule.hasTag("Accepted") && molecule.hasTag("FRET")) {
  			global_E += mol_E
  			global_S += mol_S
  			global_count += mol_obs
		}

		molecule.setParameter(E_name, mol_E / mol_obs)
		molecule.setParameter(S_name, mol_S / mol_obs)
  	}
  	builder.addParameter(E_name, global_E / global_count)
  	builder.addParameter(S_name, global_S / global_count)
}

def beta_gamma_calculation() {
	//initiate the regression
	def regression = new SimpleRegression(true)
	//add data to the model
	archive.molecules().filter{ m -> m.hasTag("Accepted") && m.hasTag("FRET")}.forEach{ molecule ->
		double iiiSapp = 0
		double iiiEapp = 0
		int count = 0
		molecule.getTable().rows().filter{ row -> !Double.isNaN(row.getValue("iiiSapp")) && !Double.isNaN(row.getValue("iiiEapp"))}.forEach{ row ->
			iiiSapp += row.getValue("iiiSapp")
			iiiEapp += row.getValue("iiiEapp")
			count++
		}
		iiiSapp = iiiSapp/count
		iiiEapp = iiiEapp/count
		regression.addData(iiiEapp,1/iiiSapp)
	}

	//save the regression outputs
	double b = regression.getSlope()
	double a = regression.getIntercept()
	double beta = a + b - 1
	double gamma = (a - 1) / (a + b - 1)

	archive.metadata().forEach{metadata ->
		metadata.setParameter("beta", beta)
		metadata.setParameter("gamma", gamma)
	}
	return [beta, gamma]
}
