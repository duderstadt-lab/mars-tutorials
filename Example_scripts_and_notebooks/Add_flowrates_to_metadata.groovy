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

	//setField(String key, double value) There is also
	//setStringField(String key, String value)...

     //FYI - Same thing
     //metadata.getPlane(0, 0, 0, t).setField(sliceToFlowMap.get(t))
     //This getPlane method directly on the MarsMetadata object takes the image index as the first input.
     //Inside it does the same as above... Maybe I should remove it?..


archive.putMetadata(metadata)
