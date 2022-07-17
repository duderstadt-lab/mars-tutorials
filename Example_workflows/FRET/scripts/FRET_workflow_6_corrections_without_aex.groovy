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

//This script was written for archives in normal memory.
//Changes are required to support virtual archives.

//This script accompanies the 'FRET dataset analysis using Mars' example pipeline as described on the mars docs.
//https://duderstadt-lab.github.io/mars-docs/examples/No_aex_FRET/

#@ MoleculeArchive archive
#@ String (label="Aem|Dex (format: channel_region)", value="532_Red") aemdexName
#@ String (label="Dem|Dex (format: channel_region)", value="532_Green") demdexName
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
import org.apache.commons.math3.fitting.GaussianCurveFitter
import org.apache.commons.math3.fitting.WeightedObservedPoints
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation

//Check that the tables of all molecule records have the columns specified
boolean foundBadRecord = false
archive.molecules().filter{molecule -> !molecule.getTable().hasColumn(aemdexName) || !molecule.getTable().hasColumn(demdexName)}\
	.findFirst().ifPresent{molecule ->
		uiService.showDialog("The molecule record " + molecule.getUID() + " is missing the Aem|Aex or Dem|Dex column specified. Aborting.\n", DialogPrompt.MessageType.ERROR_MESSAGE);
		foundBadRecord = true
}
if (foundBadRecord) return

builder = new LogBuilder()
int DO_tag_count = 0
int FRET_tag_count = 0
archive.molecules().filter{ molecule -> molecule.hasTag("Accepted")}.forEach{ molecule ->
	if (molecule.hasTag("DO"))
		DO_tag_count++
	if (molecule.hasTag("FRET"))
		FRET_tag_count++
}

//Build log message
String log = LogBuilder.buildTitleBlock("FRET workflow 6 corrections without aex")
builder.addParameter("Workflow version", "0.1")
builder.addParameter("Accepted DO", DO_tag_count)
builder.addParameter("Accepted FRET", FRET_tag_count)
builder.addParameter("Aem|Dex", aemdexName)
builder.addParameter("Dem|Dex", demdexName)

//Provide live feedback on processing in the window
if (!headless) archive.getWindow().logln("Running FRET workflow 6 corrections without aex verion 0.1")
if (!headless) archive.getWindow().logln("Found " + DO_tag_count + " accepted DO molecules.")
if (!headless) archive.getWindow().logln("Found " + FRET_tag_count + " accepted FRET molecules.")

//Background correction for each channel creating corrected columns iiIaemaex, iiIdemdex & iiIaemdex
if (!headless) archive.getWindow().logln("Running background correction and creating corrected columns iiIdemdex & iiIaemdex.")
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
	if (molecule.hasPosition("Acceptor_Bleach"))
		acceptorBleachT = molecule.getPosition("Acceptor_Bleach").getPosition()

	//BG correction for Iaemdex
	double afterBothBleach = (donorBleachT >= acceptorBleachT) ? donorBleachT : acceptorBleachT
	double background = table.mean(aemdexName, "T", afterBothBleach, table.getValue("T", table.getRowCount() - 1))
	table.rows().forEach{ row -> row.setValue("iiIaemdex", row.getValue(aemdexName) - background)}
}

// Calculation of iiEapp (E_PR)
if (!headless) archive.getWindow().logln("Calculating iiEapp and alpha.")
Calc_E_S("iiEapp", "iiIaemdex", "iiIdemdex")

// Calculation of the alpha factor
double alpha = alpha_calculation()
builder.addParameter("alpha", alpha)

// Calculation of FAD ('iiiIaemdex') values corrected for alpha
if (!headless) archive.getWindow().logln("Calculating FAD ('iiiIaemdex') values corrected for alpha.")
archive.parallelMolecules().forEach{molecule ->
    molecule.getTable().rows().forEach{ row ->
    	double Iaemdex = row.getValue("iiIaemdex")
		double Idemdex = row.getValue("iiIdemdex")
		double FAD = Iaemdex - alpha * Idemdex
		row.setValue("FAD", FAD)
    }
}

// Calculation of the beta and gamma factors
if (!headless) archive.getWindow().logln("Calculating gamma")

double gamma = gamma_calculation()
builder.addParameter("gamma", gamma)

// Calculation of FDD ("iiiIdemdex")
if (!headless) archive.getWindow().logln("Calculating FDD ('iiiIdemdex')")
archive.parallelMolecules().forEach{molecule ->
    molecule.getTable().rows().forEach{ row ->
    	double Idemdex = row.getValue("iiIdemdex")
		  double FDD = gamma * Idemdex
		  row.setValue("FDD", FDD)
    }
}

//6.F Calculation of the fully corrected S and E values
if (!headless) archive.getWindow().logln("Calculating fully corrected E.")
Calc_E_S("E", "FAD", "FDD")

//Calculating SUM_Dex (FAD+FDD) and validation parameters
if (!headless) archive.getWindow().logln("Calculating SUM_Dex (FAD+FDD) and validation parameters")
archive.parallelMolecules().forEach{molecule ->
	MarsTable table = molecule.getTable()
	for (int row=0; row < table.getRowCount(); row++) {
		double FDD = table.getValue("FDD", row)
		double FAD = table.getValue("FAD", row)
		table.setValue("SUM_Dex",row, FAD+FDD)
	}

	double lastT = table.getValue("T", table.getRowCount()-1)
	if (molecule.hasTag("FRET")) {
		double acceptorBleachT = molecule.getPosition("Acceptor_Bleach").getPosition()
		double donorBleachT = molecule.getPosition("Donor_Bleach").getPosition()
		double fretEndT = getTendFRET(molecule)
		double lastBleachT = getLastBleachT(molecule)

		double mean_SUM_Dex_FRET = table.mean("SUM_Dex", "T", 0, fretEndT)
		molecule.setParameter("SUM_Dex_FRET_Coefficient_of_Variation", table.std("SUM_Dex", "T", 0, fretEndT)/mean_SUM_Dex_FRET)

		if (acceptorBleachT < donorBleachT) {
			double mean_FDD_Donor_Recovery = table.mean("FDD", "T", acceptorBleachT, donorBleachT)
			molecule.setParameter("FDD_Donor_Recovery_Coefficient_of_Variation", table.std("FDD", "T", acceptorBleachT, donorBleachT)/mean_FDD_Donor_Recovery)
		}

		//Add E between bleaching events
		double mol_E = 0
		int mol_obs = 0
		for (int t=fretEndT; t <= lastBleachT; t++) {
  		double Iaemdex_val = table.getValue("FAD",t)
  		double Idemdex_val = table.getValue("FDD",t)
  		double E = Iaemdex_val / (Iaemdex_val + Idemdex_val)
  		if (!Double.isNaN(E)) {
    		mol_E += E
				mol_obs++
  		}
		}
		molecule.setParameter("E_Between_Bleaches", mol_E / mol_obs)

		double[] fret_FDD = table.getColumnAsDoublesNoNaNs("FDD", "T", 0, fretEndT)
		double[] fret_FAD = table.getColumnAsDoublesNoNaNs("FAD", "T", 0, fretEndT)
		try {
			molecule.setParameter("FRET_Pearsons_Correlation", new PearsonsCorrelation().correlation(fret_FDD, fret_FAD))
		} catch(Exception e) {
			molecule.setParameter("FRET_Pearsons_Correlation", Double.NaN)
		}
	}

	if (molecule.hasTag("DO")) {
		double donorBleachT = molecule.getPosition("Donor_Bleach").getPosition()
		double mean_FDD = table.mean("FDD", "T", 0, donorBleachT)
		molecule.setParameter("FDD_Coefficient_of_Variation", table.std("FDD", "T", 0, donorBleachT)/mean_FDD)
	}
}

if (!headless) archive.getWindow().logln("Building final log.")
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

def getLastBleachT(def molecule) {
	double donorBleachT = 0
	double acceptorBleachT = 0
	if (molecule.hasPosition("Donor_Bleach"))
		donorBleachT = molecule.getPosition("Donor_Bleach").getPosition()
	if (molecule.hasPosition("Acceptor_Bleach"))
		acceptorBleachT = molecule.getPosition("Acceptor_Bleach").getPosition()

	return (donorBleachT >= acceptorBleachT) ? donorBleachT : acceptorBleachT
}

def Calc_E_S(String E_name, String Iaemdex, String Idemdex) {
	double global_E = 0
	int global_count = 0
	archive.molecules().forEach{molecule ->
		MarsTable table = molecule.getTable()
  		double endT = 0
  		if (molecule.hasTag("DO"))
  			endT = molecule.getPosition("Donor_Bleach").getPosition()
		else if (molecule.hasTag("FRET"))
			endT = getTendFRET(molecule)

  		double mol_E = 0
		  int mol_obs = 0
  		for (int i=0; i < endT; i++) {
    		double Iaemdex_val = table.getValue(Iaemdex,i)
    		double Idemdex_val = table.getValue(Idemdex,i)
    		double E = Iaemdex_val / (Iaemdex_val + Idemdex_val)
    		if (!Double.isNaN(E)) {
	    		table.setValue(E_name,i,E)
	    		mol_E += E
					mol_obs++
    		}
  		}

		if (molecule.hasTag("Accepted") && molecule.hasTag("FRET")) {
  			global_E += mol_E / mol_obs
  			global_count++
		}

		molecule.setParameter(E_name, mol_E / mol_obs)
  }
  builder.addParameter(E_name, global_E / global_count)
}

def alpha_calculation() {
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

	//Negative values are not possible. If negative we set to zero
	double alpha = E_stat/(1-E_stat)
	if (alpha < 0) alpha = 0

	archive.metadata().forEach{metadata ->
		metadata.setParameter("iiEappDO_mean", E_stat)
		metadata.setParameter("alpha", alpha)
	}

	return alpha
}

//Only FRET traces where the donor survives more than 25 time points after acceptor bleaching
//are used for the gamma calculation.
//Only gamma values from individual molecules between 0 and 5 are used to determine the global gamma value
def gamma_calculation() {
	//No great histogram function in groovy so we do it outselves :(
	//We calculate the sqrt of the observations for the bin count in the histogram
	int bins = (int)Math.sqrt(archive.molecules().filter{ m -> m.hasTag("Accepted") && m.hasTag("FRET")}\
	.filter{ m -> m.getPosition("Donor_Bleach").getPosition() > m.getPosition("Acceptor_Bleach").getPosition()}\
	.filter{ m -> m.getPosition("Donor_Bleach").getPosition() - m.getPosition("Acceptor_Bleach").getPosition() > 25}.count())
	double minX = 0
	double maxX = 5
	double binWidth = (maxX - minX) / bins
	double[] yvalues = new double[bins]
	double[] xvalues = new double[bins]

	for (int bin = 0; bin < bins; bin++) {
		yvalues[bin] = 0
		xvalues[bin] = minX + (0.5 + bin) * binWidth
	}

	//add data to the model
	archive.molecules().filter{ m -> m.hasTag("Accepted") && m.hasTag("FRET")}\
	.filter{ m -> m.getPosition("Donor_Bleach").getPosition() > m.getPosition("Acceptor_Bleach").getPosition()}\
	.filter{ m -> m.getPosition("Donor_Bleach").getPosition() - m.getPosition("Acceptor_Bleach").getPosition() > 25}\
	.forEach{ molecule ->
		MarsTable table = molecule.getTable()

		double iAPre =  table.mean(aemdexName, "T", 0, molecule.getPosition("Acceptor_Bleach").getPosition())
		double iAPost = table.mean(aemdexName, "T", molecule.getPosition("Acceptor_Bleach").getPosition(), molecule.getPosition("Donor_Bleach").getPosition())

		double iDPre =  table.mean(demdexName, "T", 0, molecule.getPosition("Acceptor_Bleach").getPosition())
		double iDPost = table.mean(demdexName, "T", molecule.getPosition("Acceptor_Bleach").getPosition(), molecule.getPosition("Donor_Bleach").getPosition())

		double gamma = (iAPre - iAPost)/(iDPost - iDPre)
		if (!Double.isNaN(gamma) && gamma > minX && gamma < maxX) {
			for (int bin = 0; bin < bins; bin++) {
				if (gamma >= minX + bin * binWidth && gamma < minX + (bin + 1) *
					binWidth)
				{
					yvalues[bin]++
					break
				}
			}
		}
		molecule.setParameter("gamma", gamma)
	}

	WeightedObservedPoints obs = new WeightedObservedPoints()
	for (int bin = 0; bin < bins; bin++)
		obs.add(xvalues[bin],  yvalues[bin])

	double[] parameters = GaussianCurveFitter.create().fit(obs.toList())
	double globalGamma = parameters[1]
	archive.metadata().forEach{metadata -> metadata.setParameter("gamma", globalGamma)}
	return globalGamma
}
