//Adds a time column to the Metadata by converting the slice number together
//with a time/slice input to a time measurement in seconds.

#@ MoleculeArchive archive
#@ double timePerSlice

import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import org.scijava.table.*

archive.lock()

//Retrieves the first image metadata item
MarsMetadata metadata = archive.getMetadata(0)
MarsTable metaDataTable = metadata.getDataTable()

DoubleColumn timeCol = new DoubleColumn("Time (s)")

for (int row=0; row<metaDataTable.getRowCount(); row++) {
   int slice = metaDataTable.getValue("slice", row);
   timeCol.add((slice-1)*timePerSlice)
}

while (metaDataTable.hasColumn("Time (s)")) {
    metaDataTable.removeColumn("Time (s)");
}

metaDataTable.add(timeCol)

archive.putMetadata(metadata)

archive.unlock()
