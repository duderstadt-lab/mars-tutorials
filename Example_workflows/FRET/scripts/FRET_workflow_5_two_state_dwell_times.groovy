#@ String (label="Efficiency (E)", value="E") efficiencyColumn
#@ String (label="Time (s)", value="532_Green_Time_(s)") timeColumn
#@ Double (label="Efficiency threshold (e.g. 0.5)", value = 0.4) threshold
#@ MoleculeArchive archive

import de.mpg.biochem.mars.table.*

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
