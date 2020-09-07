//Adds flowrates from metadata to molecule records
//Run 'Add flowrate (ul/min) to the metadata' first to add flowrates to the metadata of the Archive

#@ MoleculeArchive archive

import de.mpg.biochem.mars.metadata.*
import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import org.scijava.table.*


MarsMetadata metadata = archive.getMetadata(0)

archive.molecules().forEach{ molecule ->
  MarsTable table = molecule.getDataTable()

  if (!table.hasColumn("Flowrate (ul/min)"))
    table.appendColumn("Flowrate (ul/min)")

  for (int row=0;row<table.getRowCount();row++){
    int T = (int) table.getValue("T", row)
    table.set("Flowrate (ul/min)", row,metadata.getImage(0).getPlane(0, 0, T).getField("Flowrate (ul/min)") )
  }
}
