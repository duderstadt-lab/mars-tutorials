#@ MoleculeArchive archive
#@ double timePerSlice

import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import org.scijava.table.*

archive.lock()

//Retrieves the first image metadata item, can instead use metaUID if needed
MarsImageMetadata metaData = archive.getImageMetadata(0)
MarsTable metaDataTable = metaData.getDataTable()

DoubleColumn timeCol = new DoubleColumn("Time (s)")

for (int row=0; row<metaDataTable.getRowCount(); row++) {
   int slice = metaDataTable.getValue("slice", row);
   timeCol.add((slice-1)*timePerSlice)
}

while (metaDataTable.hasColumn("Time (s)")) {
    metaDataTable.removeColumn("Time (s)");
}

metaDataTable.add(timeCol)

archive.putImageMetadata(metaData)

archive.unlock()
