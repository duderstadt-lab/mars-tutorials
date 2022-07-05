/*******************************************************************************
 * Copyright (C) 2022, Duderstadt Lab
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
//written by: Nadia M. Huisjes, MSc. & Karl E. Duderstadt

//This script accompanies the 'FRET dataset analysis using Mars' example pipeline as described on the mars docs.
//https://duderstadt-lab.github.io/mars-docs/examples/FRET

#@ String (label="Aem|Aex (format: channel_region)", value="637_Red_Profile_Corrected") aemaexName
#@ String (label="Aem|Dex (format: channel_region)", value="532_Red") aemdexName
#@ String (label="Dem|Dex (format: channel_region)", value="532_Green") demdexName
#@ String (label="Gamma model", choices={"static molecules", "dynamic molecules"}, value="dynamic molecules", style="radioButtonVertical") gammaModel
#@ Double (label="Dynamic efficiency threshold", value = 0.4) gammaEfficiencyThreshold
#@ MoleculeArchive archive
#@ UIService uiService

headless = (archive.getWindow() == null) ? true : false

//Import packages
import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import de.mpg.biochem.mars.util.*
import org.scijava.table.*
import groovy.lang.*
import de.mpg.biochem.mars.kcp.commands.*
import org.scijava.ui.DialogPrompt
import org.apache.commons.math3.stat.regression.SimpleRegression

//Check that the tables of all molecule records have the columns specified
boolean foundBadRecord = false
archive.molecules().filter{molecule -> !molecule.getTable().hasColumn(aemaexName) || !molecule.getTable().hasColumn(aemdexName) || !molecule.getTable().hasColumn(demdexName)}\
	.findFirst().ifPresent{molecule ->
		uiService.showDialog("The molecule record " + molecule.getUID() + " is missing the Aem|Aex, Aem|Dex or Dem|Dex column specified. Aborting.\n", DialogPrompt.MessageType.ERROR_MESSAGE);
		foundBadRecord = true
}
if (foundBadRecord) return

builder = new LogBuilder()
int DO_tag_count = 0
int AO_tag_count = 0
int FRET_tag_count = 0
archive.molecules().filter{ molecule -> molecule.hasTag("Accepted")}.forEach{ molecule ->
	if (molecule.hasTag("DO"))
		DO_tag_count++
	if (molecule.hasTag("AO"))
		AO_tag_count++
	if (molecule.hasTag("FRET"))
		FRET_tag_count++
}

//Build log message
String log = LogBuilder.buildTitleBlock("FRET workflow 4 alex corrections")
builder.addParameter("Workflow version", "0.1")
builder.addParameter("Accepted DO", DO_tag_count)
builder.addParameter("Accepted AO", AO_tag_count)
builder.addParameter("Accepted FRET", FRET_tag_count)
builder.addParameter("Aem|Aex", aemaexName)
builder.addParameter("Aem|Dex", aemdexName)
builder.addParameter("Dem|Dex", demdexName)
builder.addParameter("Gamma model", gammaModel)
builder.addParameter("Dynamic efficiency threshold", gammaEfficiencyThreshold)

//Provide live feedback on processing in the window
if (!headless) archive.getWindow().logln("Running FRET alex corrections workflow version 0.1")
if (!headless) archive.getWindow().logln("Found " + DO_tag_count + " accepted DO molecules.")
if (!headless) archive.getWindow().logln("Found " + AO_tag_count + " accepted AO molecules.")
if (!headless) archive.getWindow().logln("Found " + FRET_tag_count + " accepted FRET molecules.")

//Background correction for each channel creating corrected columns iiIaemaex, iiIdemdex & iiIaemdex
if (!headless) archive.getWindow().logln("Running background correction for each channel and creating corrected columns iiIaemaex, iiIdemdex & iiIaemdex.")
archive.parallelMolecules().forEach{molecule ->
	MarsTable table = molecule.getTable()

	double donorBleachT = 0
	double acceptorBleachT = 0

	//BG correction for Idemdex
	if (molecule.hasPosition("Donor_Bleach")) {
		donorBleachT = molecule.getPosition("Donor_Bleach").getPosition()
		double donorBackground = table.mean(demdexName, "T", donorBleachT, table.getValue("T", table.getRowCount() - 1))
		table.rows().forEach{ row -> row.setValue("iiIdemdex", row.getValue(demdexName) - donorBackground)}
	} else {
		double donorBackground = table.mean(demdexName, "T", 0, table.getValue("T", table.getRowCount() - 1))
		table.rows().forEach{ row -> row.setValue("iiIdemdex", row.getValue(demdexName) - donorBackground)}
	}

	//BG correction for Iaemaex
	if (molecule.hasPosition("Acceptor_Bleach")) {
		acceptorBleachT = molecule.getPosition("Acceptor_Bleach").getPosition()
		double acceptorBackground = table.mean(aemaexName, "T", acceptorBleachT, table.getValue("T", table.getRowCount() - 1))
		table.rows().forEach{ row -> row.setValue("iiIaemaex", row.getValue(aemaexName) - acceptorBackground)}
	} else {
		double acceptorBackground = table.mean(aemaexName, "T", 0, table.getValue("T", table.getRowCount() - 1))
		table.rows().forEach{ row -> row.setValue("iiIaemaex", row.getValue(aemaexName) - acceptorBackground)}
	}

	//BG correction for Iaemdex
	double afterBothBleach = (donorBleachT >= acceptorBleachT) ? donorBleachT : acceptorBleachT
	double background = table.mean(aemdexName, "T", afterBothBleach, table.getValue("T", table.getRowCount() - 1))
	table.rows().forEach{ row -> row.setValue("iiIaemdex", row.getValue(aemdexName) - background)}
}

//5.C Calculation of iiEapp and iiSapp values
if (!headless) archive.getWindow().logln("Calculating iiEapp and iiSapp and alpha and delta.")
Calc_E_S("iiEapp", "iiSapp", "iiIaemdex", "iiIdemdex", "iiIaemaex")

//6. Data Corrections alpha, beta, gamma, delta
//6.A Calculation of the alpha and delta factors
alpha_delta_list = alpha_delta_calculation()
double alpha = alpha_delta_list[0]
double delta = alpha_delta_list[1]

builder.addParameter("alpha", alpha)
builder.addParameter("delta", delta)

//6.B Calculation of FAD ('iiiIaemdex') values corrected for alpha and delta
if (!headless) archive.getWindow().logln("Calculating FAD ('iiiIaemdex') values corrected for alpha and delta.")
archive.parallelMolecules().forEach{molecule ->
    molecule.getTable().rows().forEach{ row ->
		double Iaemaex = row.getValue("iiIaemaex")
		double Iaemdex = row.getValue("iiIaemdex")
		double Idemdex = row.getValue("iiIdemdex")
		double FAD = Iaemdex - alpha * Idemdex - delta * Iaemaex
		row.setValue("FAD", FAD)
    }
}

//6.C Calculation of iiiSapp and iiiEapp (corrected for alpha and delta)
if (!headless) archive.getWindow().logln("Calculating iiiEapp and iiiSapp corrected for alpha and delta.")
Calc_E_S("iiiEapp", "iiiSapp", "FAD", "iiIdemdex", "iiIaemaex")

//6.D Calculation of the beta and gamma factors
if (!headless) archive.getWindow().logln("Calculating beta and gamma for " + gammaModel)
def beta_gamma = (gammaModel.equals("dynamic molecules")) ? dynamic_molecules_beta_gamma_calculation() : static_molecules_beta_gamma_calculation()
double beta = beta_gamma[0]
double gamma = beta_gamma[1]

builder.addParameter("beta", beta)
builder.addParameter("gamma", gamma)

//6.E Calculation of FAA ("iiiIaemaex") and FDD ("iiiIdemdex")
if (!headless) archive.getWindow().logln("Calculating FAA ('iiiIaemaex') and FDD ('iiiIdemdex')")
archive.parallelMolecules().forEach{molecule ->
    molecule.getTable().rows().forEach{ row ->
		double Iaemaex = row.getValue("iiIaemaex")
		double Idemdex = row.getValue("iiIdemdex")
		double FDD = gamma * Idemdex
		double FAA = Iaemaex / beta
		row.setValue("FDD",FDD)
		row.setValue("FAA",FAA)
    }
}

//Calculating SUM_Dex (FAD+FDD) and SUM_signal (FAA+FAD+FDD)
if (!headless) archive.getWindow().logln("Calculating SUM_Dex (FAD+FDD) and SUM_signal (FAA+FAD+FDD)")
archive.parallelMolecules().forEach{molecule -> 
	MarsTable table = molecule.getTable()
	for (int row=0; row < table.getRowCount(); row++) {
		double FDD = table.getValue("FDD", row)
		double FAA = table.getValue("FAA", row)
		double FAD = table.getValue("FAD", row)
		table.setValue("SUM_Dex",row, FAD+FDD)
		table.setValue("SUM_signal",row, FAA+FAD+FDD)
	}
}

//6.F Calculation of the fully corrected S and E values
if (!headless) archive.getWindow().logln("Calculating fully corrected E and S. Building final log.")
Calc_E_S("E", "S", "FAD", "FDD", "FAA")

log += builder.buildParameterList()
log += "\n" + LogBuilder.endBlock(true) + "\n"
archive.logln(log)
if (!headless) archive.getWindow().logln("Done!")

//Utility methods
def getTendFRET(def molecule) {
	double donorBleachT = 0
	double acceptorBleachT = 0
	if (molecule.hasPosition("Donor_Bleach"))
		donorBleachT = molecule.getPosition("Donor_Bleach").getPosition()
	if (molecule.hasPosition("Acceptor_Bleach"))
		acceptorBleachT = molecule.getPosition("Acceptor_Bleach").getPosition()

	return (donorBleachT <= acceptorBleachT) ? donorBleachT : acceptorBleachT
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
  		for (int i=0; i < endT; i++) {
    		double Iaemaex_val = table.getValue(Iaemaex,i)
    		double Iaemdex_val = table.getValue(Iaemdex,i)
    		double Idemdex_val = table.getValue(Idemdex,i)
  			double S = (Iaemdex_val + Idemdex_val) / (Iaemdex_val + Idemdex_val + Iaemaex_val)
    		double E = Iaemdex_val / (Iaemdex_val + Idemdex_val)
    		if (!Double.isNaN(S) && !Double.isNaN(E)) {
	    		table.setValue(E_name,i,E)
	    		table.setValue(S_name,i,S)
	    		mol_E += E
	    		mol_S += S
					mol_obs++
    		}
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

def alpha_delta_calculation() {
	double observations = 0
	double E_stat = 0
	archive.molecules().filter{ m -> m.hasTag("DO")}.forEach{ molecule ->
		double mol_E_stat = 0
		int mol_obs = 0
		molecule.getTable().rows().filter{ r -> !Double.isNaN(r.getValue("iiEapp"))}.forEach{ row ->
			mol_E_stat += row.getValue("iiEapp")
			mol_obs++
		}
		mol_E_stat = mol_E_stat / mol_obs

		if (molecule.hasTag("Accepted")) {
			E_stat += mol_E_stat
			observations++
		}
		molecule.setParameter("alpha", mol_E_stat/(1-mol_E_stat))
	}
	E_stat = E_stat / observations

	observations = 0
	double S_stat = 0
	archive.molecules().filter{ m -> m.hasTag("AO")}.forEach{ molecule ->
		double mol_S_stat = 0
		double mol_obs = 0
		molecule.getTable().rows().filter{ r -> !Double.isNaN(r.getValue("iiSapp"))}.forEach{ row ->
			mol_S_stat += row.getValue("iiSapp")
			mol_obs++
		}
		mol_S_stat = mol_S_stat / mol_obs

		if (molecule.hasTag("Accepted")) {
			S_stat += mol_S_stat
			observations++
		}
		molecule.setParameter("delta", mol_S_stat/(1-mol_S_stat))
	}
	S_stat = S_stat / observations

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

def static_molecules_beta_gamma_calculation() {
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

def dynamic_molecules_beta_gamma_calculation() {
	//initiate the regression
	def regression = new SimpleRegression(true)
	//add data to the model
	archive.molecules().filter{ m -> m.hasTag("Accepted") && m.hasTag("FRET")}.forEach{ molecule ->
		double iiiSappHigh = 0
		double iiiEappHigh = 0
		double countHigh = 0
		double iiiSappLow = 0
		double iiiEappLow = 0
		double countLow = 0
		molecule.getTable().rows().filter{ row -> !Double.isNaN(row.getValue("iiiSapp")) && !Double.isNaN(row.getValue("iiiEapp"))}.forEach{ row ->
			if (row.getValue("iiiEapp") > gammaEfficiencyThreshold) {
				iiiSappHigh += row.getValue("iiiSapp")
				iiiEappHigh += row.getValue("iiiEapp")
				countHigh++
			} else {
				iiiSappLow += row.getValue("iiiSapp")
				iiiEappLow += row.getValue("iiiEapp")
				countLow++
			}
		}
		iiiSappHigh = iiiSappHigh/countHigh
		iiiEappHigh = iiiEappHigh/countHigh
		if (!Double.isNaN(iiiSappHigh) && !Double.isNaN(iiiEappHigh)) regression.addData(iiiEappHigh,1/iiiSappHigh)
		iiiSappLow = iiiSappLow/countLow
		iiiEappLow = iiiEappLow/countLow
		if (!Double.isNaN(iiiSappLow) && !Double.isNaN(iiiEappLow)) regression.addData(iiiEappLow,1/iiiSappLow)
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
