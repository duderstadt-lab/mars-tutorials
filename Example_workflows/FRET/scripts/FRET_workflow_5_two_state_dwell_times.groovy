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

//This script accompanies the 'FRET dataset analysis using Mars' example pipeline as described on the mars docs.
//https://duderstadt-lab.github.io/mars-docs/examples/FRET_dynamic

#@ String (label="Efficiency (E)", value="E") efficiencyColumn
#@ String (label="Time (s)", value="532_Green_Time_(s)") timeColumn
#@ Double (label="Efficiency threshold (e.g. 0.5)", value = 0.4) threshold
#@ MoleculeArchive archive

import de.mpg.biochem.mars.table.*
import de.mpg.biochem.mars.util.*

//Build log message
builder = new LogBuilder()
String log = LogBuilder.buildTitleBlock("FRET two state dwell times workflow 5")
builder.addParameter("Workflow version", "0.1")
builder.addParameter("Efficiency (E)", efficiencyColumn)
builder.addParameter("Time (s)", timeColumn)
builder.addParameter("Efficiency threshold (e.g. 0.5)", threshold)
log += builder.buildParameterList()

//First calculate mean populations
double highE = 0
int highEcount = 0
double lowE = 0
int lowEcount = 0

archive.molecules().filter{ m -> m.hasTag("FRET") && m.hasTag("Accepted")}.forEach{ molecule ->
	MarsTable table = molecule.getTable()
	table.rows().filter{ row -> !Double.isNaN(row.getValue(efficiencyColumn))}.forEach{ row ->
		double e = row.getValue(efficiencyColumn)
		if (e > threshold) {
			highE += e
			highEcount++
		} else {
			lowE += e
			lowEcount++
		}
	}
}

double meanHighE = highE / highEcount
double meanLowE = lowE / lowEcount

archive.molecules().filter{ m -> m.hasTag("FRET") && m.hasTag("Accepted")}.forEach{ molecule ->
	MarsTable segmentsTable = new MarsTable(efficiencyColumn + " vs " + timeColumn, "X1", "Y1", "X2", "Y2", "A", "Sigma_A", "B", "Sigma_B")

	MarsTable table = molecule.getTable()
	double prevE = table.getValue(efficiencyColumn, 0)
	double x1 = table.getValue(timeColumn,0)
	for (int row = 1; row < table.getRowCount(); row++) {
		double e = table.getValue(efficiencyColumn, row)
		if (Double.isNaN(e)) {
			addSegmentsTableRow(segmentsTable, x1, table.getValue(timeColumn, row), (prevE > threshold) ? meanHighE : meanLowE)
			break
		}

		//Check for state switching
		if (prevE > threshold && e < threshold) {
			//high to low switch
			addSegmentsTableRow(segmentsTable, x1, table.getValue(timeColumn, row), meanHighE)
			x1 = table.getValue(timeColumn, row)
		} else if (prevE < threshold && e > threshold) {
			//low to high switch
			addSegmentsTableRow(segmentsTable, x1, table.getValue(timeColumn, row), meanLowE)
			x1 = table.getValue(timeColumn, row)
		}
		prevE = e
	}
	molecule.putSegmentsTable(timeColumn, efficiencyColumn, segmentsTable)
}

def addSegmentsTableRow(def sTable, def rx1, def rx2, def re) {
	int r = sTable.getRowCount()
	sTable.appendRow()
	sTable.setValue("X1", r, rx1)
	sTable.setValue("Y1", r, re)
	sTable.setValue("X2", r, rx2)
	sTable.setValue("Y2", r, re)
	sTable.setValue("A", r, re)
	sTable.setValue("Sigma_A", r, 0)
	sTable.setValue("B", r, 0)
	sTable.setValue("Sigma_B", r, 0)
}
log += "\n" + LogBuilder.endBlock(true) + "\n"
archive.logln(log)
