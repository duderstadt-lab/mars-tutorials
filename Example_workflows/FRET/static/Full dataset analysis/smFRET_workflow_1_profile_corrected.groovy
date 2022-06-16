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
import groovy.lang.*

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

String log = LogBuilder.buildTitleBlock("smFRET workflow 1")
builder.addParameter("Workflow version", "0.2")
builder.addParameter("Accepted DO", DO_tag_count)
builder.addParameter("Accepted AO", AO_tag_count)
builder.addParameter("Accepted FRET", FRET_tag_count)

//Background correction for each channel creating corrected columns iiIaemaex, iiIdemdex & iiIaemdex
archive.molecules().forEach{molecule ->
	MarsTable table = molecule.getTable()

	double donorBleachT = 0
	double acceptorBleachT = 0

	//BG correction for Idemdex
	if (molecule.hasPosition("Donor_Bleach")) {
		donorBleachT = molecule.getPosition("Donor_Bleach").getPosition()
		double donorBackground = table.mean("1_Green", "T", donorBleachT, table.getValue("T", table.getRowCount() - 1))
		table.rows().forEach{ row -> row.setValue("iiIdemdex", row.getValue("1_Green") - donorBackground)}
	} else {
		double donorBackground = table.mean("1_Green", "T", 0, table.getValue("T", table.getRowCount() - 1))
		table.rows().forEach{ row -> row.setValue("iiIdemdex", row.getValue("1_Green") - donorBackground)}
	}

	//BG correction for Iaemaex
	if (molecule.hasPosition("Acceptor_Bleach")) {
		acceptorBleachT = molecule.getPosition("Acceptor_Bleach").getPosition()
		double acceptorBackground = table.mean("0_Profile_Corrected", "T", acceptorBleachT, table.getValue("T", table.getRowCount() - 1))
		table.rows().forEach{ row -> row.setValue("iiIaemaex", row.getValue("0_Profile_Corrected") - acceptorBackground)}
	} else {
		double acceptorBackground = table.mean("0_Profile_Corrected", "T", 0, table.getValue("T", table.getRowCount() - 1))
		table.rows().forEach{ row -> row.setValue("iiIaemaex", row.getValue("0_Profile_Corrected") - acceptorBackground)}
	}

	//BG correction for Iaemdex
	double afterBothBleach = (donorBleachT >= acceptorBleachT) ? donorBleachT : acceptorBleachT
	double background = table.mean("1_Red", "T", afterBothBleach, table.getValue("T", table.getRowCount() - 1))
	table.rows().forEach{ row -> row.setValue("iiIaemdex", row.getValue("1_Red") - background)}
}

//5.C Calculation of iiEapp and iiSapp values
Calc_E_S("iiEapp", "iiSapp", "iiIaemdex", "iiIdemdex", "iiIaemaex")

//6. Data Corrections alpha, beta, gamma, delta
//6.A Calculation of the alpha and delta factors
alpha_delta_list = alpha_delta_calculation()
double alpha = alpha_delta_list[0]
double delta = alpha_delta_list[1]

//6.B Calculation of FAD ('iiiIaemdex') values corrected for alpha and delta
archive.molecules().forEach{molecule ->
    MarsTable table = molecule.getTable()
    double endT = 0
	if (molecule.hasTag("DO"))
		endT = molecule.getPosition("Donor_Bleach").getPosition()
	else if (molecule.hasTag("AO"))
		endT = molecule.getPosition("Acceptor_Bleach").getPosition()
	else if (molecule.hasTag("FRET"))
		endT = getTendFRET(molecule)

    double len = table.getRowCount()
	for (i=0; i<endT; i++){
		double Iaemaex = table.getValue("iiIaemaex",i)
		double Iaemdex = table.getValue("iiIaemdex",i)
		double Idemdex = table.getValue("iiIdemdex",i)
		double FAD = Iaemdex - alpha * Idemdex - delta * Iaemaex
		table.setValue("FAD",i,FAD)
  	}
  	for (i=endT; i<len; i++){
  		table.setValue("FAD",(int)i,"NaN")
    }
}

//6.C Calculation of iiiSapp and iiiEapp (corrected for alpha and delta)
Calc_E_S("iiiEapp", "iiiSapp", "FAD", "iiIdemdex", "iiIaemaex")

builder.addParameter("alpha", alpha)
builder.addParameter("delta", delta)
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

def alpha_delta_calculation() {
	int observations = 0
	double E_stat = 0
	archive.molecules().filter{ m -> m.hasTag("DO")}.forEach{ molecule ->
		double mol_E_stat = 0
		int mol_obs = 0
		molecule.getTable().rows().filter{ r -> !Double.isNaN(r.getValue("iiEapp"))}.forEach{ row ->
			mol_E_stat += row.getValue("iiEapp")
			mol_obs++
		}
		if (molecule.hasTag("Accepted")) {
			E_stat += mol_E_stat
			observations += mol_obs
		}

		mol_E_stat = mol_E_stat / mol_obs
		molecule.setParameter("alpha", mol_E_stat/(1-mol_E_stat))
	}
	E_stat = E_stat / observations

	observations = 0
	double S_stat = 0
	archive.molecules().filter{ m -> m.hasTag("AO")}.forEach{ molecule ->
		double mol_S_stat = 0
		int mol_obs = 0
		molecule.getTable().rows().filter{ r -> !Double.isNaN(r.getValue("iiSapp"))}.forEach{ row ->
			mol_S_stat += row.getValue("iiSapp")
			mol_obs++
		}
		if (molecule.hasTag("Accepted")) {
			S_stat += mol_S_stat
			observations += mol_obs
		}

		mol_S_stat = mol_S_stat / mol_obs
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