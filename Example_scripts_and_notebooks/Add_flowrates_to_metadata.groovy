//Adds a flowrate (ul/min) column to the metadata record.

#@ MarsTable sliceToFlowTable
#@ MoleculeArchive archive
#@ String metaUID

import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import org.scijava.table.*
import java.util.HashMap;

//First lets build a map from slice to flowrate
sliceToFlowMap = [:]

for (int row=0;row<sliceToFlowTable.getRowCount();row++)
   sliceToFlowMap.put(sliceToFlowTable.getValue("T", row), sliceToFlowTable.getValue("flowrate", row))
newColumn = "Flowrate (ul/min)"

MarsMetadata metaData = archive.getMetadata(0)
MarsTable metaTable = metaData.getDataTable()

if (!metaTable.hasColumn(newColumn))
  metaTable.appendColumn(newColumn)

for (int row=0;row<metaTable.getRowCount();row++) {
  if (sliceToFlowMap.containsKey(metaTable.getValue("T", row)))
      metaTable.setValue(newColumn, row, sliceToFlowMap.get(metaTable.getValue("T", row)))
}

archive.putMetadata(metaData)

//This requires the input of a MarsTable with flow rates. A dummy example is created with this script.
#@output MarsTable table

import de.mpg.biochem.mars.table.*;

//Initialize a new ResultsTable with 2 columns and no rows
table = new MarsTable(2,0)
table.setColumnHeader(0, "T")
table.setColumnHeader(1, "flowrate")


for (int row=0;row<150;row++) {
   //Since we initialised the table with no rows
   //We have to increase the table size by one row
   //before we can add more values
   table.appendRow()
   table.setValue("T", row, row)
   table.setValue("flowrate", row, row*0.25)

}
