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
//written by: Karl E. Duderstadt

//This script was written for archives in normal memory.
//Changes are required to support virtual archives.

//This script accompanies the 'FRET dataset analysis using Mars' example workflows as described on the mars docs:
//https://duderstadt-lab.github.io/mars-docs/examples/Static_FRET and
//https://duderstadt-lab.github.io/mars-docs/examples/Dynamic_FRET and
//https://duderstadt-lab.github.io/mars-docs/examples/No_aex_FRET/

#@ MoleculeArchive archive
#@ Double (label="SUM_Dex_FRET_Coefficient_of_Variation < ", style="slider", min=0, max=1, stepSize=0.01, value=0.4) SUM_Dex_FRET_Coefficient_of_Variation
#@ Double (label="SUM_signal_FRET_Coefficient_of_Variation < ", style="slider", min=0, max=1, stepSize=0.01, value=0.3) SUM_signal_FRET_Coefficient_of_Variation
#@ Double (label="FRET_Pearsons_Correlation < ", style="slider", min=-1, max=1, stepSize=0.01, value=-0.3) FRET_Pearsons_Correlation
#@ String (label="Tagging mode:", choices={"Tag valid molecules with Accepted", "Remove Accepted tag from invalid molecules"}, style="radioButtonVertical") tagMode

import de.mpg.biochem.mars.util.*

//Build log message
builder = new LogBuilder()
String log = LogBuilder.buildTitleBlock("FRET workflow 7 validation filter")
builder.addParameter("Workflow version", "0.1")
builder.addParameter("SUM_Dex_FRET_Coefficient_of_Variation < ", SUM_Dex_FRET_Coefficient_of_Variation)
builder.addParameter("SUM_signal_FRET_Coefficient_of_Variation < ", SUM_signal_FRET_Coefficient_of_Variation)
builder.addParameter("FRET_Pearsons_Correlation < ", FRET_Pearsons_Correlation)
log += builder.buildParameterList()
archive.logln(log)

archive.molecules().forEach{ molecule ->
	boolean valid = true
	if (molecule.hasParameter("SUM_Dex_FRET_Coefficient_of_Variation") && molecule.getParameter("SUM_Dex_FRET_Coefficient_of_Variation") > SUM_Dex_FRET_Coefficient_of_Variation)
		valid = false

	if (molecule.hasParameter("SUM_signal_FRET_Coefficient_of_Variation") && molecule.getParameter("SUM_signal_FRET_Coefficient_of_Variation") > SUM_signal_FRET_Coefficient_of_Variation)
		valid = false

	if (molecule.hasParameter("FRET_Pearsons_Correlation") && molecule.getParameter("FRET_Pearsons_Correlation") > FRET_Pearsons_Correlation)
		valid = false

	if (tagMode.equals("Tag valid molecules with Accepted") && valid)
		molecule.addTag("Accepted")
	else if (tagMode.equals("Remove Accepted tag from invalid molecules") && !valid)
		molecule.removeTag("Accepted")
}

archive.logln(LogBuilder.endBlock(true))
