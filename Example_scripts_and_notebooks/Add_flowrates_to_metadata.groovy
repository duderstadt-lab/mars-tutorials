//Adds a flowrate (ul/min) column to the metadata record.

#@ MarsTable tToFlowTable
#@ MoleculeArchive archive

import de.mpg.biochem.mars.metadata.*
import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import org.scijava.table.*
import java.util.HashMap;

//First lets build a map from slice to flowrate
tToFlowMap = [:]

for (int row=0;row<tToFlowTable.getRowCount();row++)
   tToFlowMap.put((int) tToFlowTable.getValue("T", row), tToFlowTable.getValue("flowrate", row))

MarsMetadata metadata = archive.getMetadata(0)

for (int indexT=0; indexT<metadata.getImage(0).getSizeT(); indexT++)
     metadata.getImage(0).getPlane(0, 0, indexT).setField("Flowrate (ul/min)", tToFlowMap.get(indexT))


archive.putMetadata(metadata)


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
