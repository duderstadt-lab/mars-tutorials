//Adds a flowrate (ul/min) column to the metadata record.

#@ MarsTable sliceToFlowTable
#@ MoleculeArchive archive
#@ String metaUID

import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import org.scijava.table.*

import java.util.HashMap;

//First lets build a map from slice to flowrate
HashMap<Double, Double> sliceToFlowMap = new HashMap<Double, Double>();
for (int row=0;row<sliceToFlowTable.getRowCount();row++) {
   sliceToFlowMap.put(sliceToFlowTable.getValue("slice", row), sliceToFlowTable.getValue("flowrate", row));
}

newColumn = "Flowrate (ul/min)"

archive.lock()

MarsImageMetadata metaData = archive.getImageMetadata(metaUID)
MarsTable metaTable = metaData.getDataTable()

if (!metaTable.hasColumn(newColumn))
  metaTable.appendColumn(newColumn)

for (int row=0;row<metaTable.getRowCount();row++)
  metaTable.setValue(newColumn, row, (double)sliceToFlowMap.get(metaTable.getValue("slice", row)))

archive.putImageMetadata(metaData);

archive.unlock()
