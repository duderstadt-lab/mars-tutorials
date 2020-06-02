//Adds flowrates from metadata to molecule records

#@ MoleculeArchive archive

import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import org.scijava.table.*

archive.lock()

MarsTable metaTable = archive.getMetadata(0).getDataTable()

sliceToFlowrate = [:]

for (int row=0;row<metaTable.getRowCount();row++) {
   sliceToFlowrate.put(metaTable.getValue("slice", row), metaTable.getValue("Flowrate (ul/min)", row))
}

println(sliceToFlowrate)

archive.getMoleculeUIDs().stream().forEach{ UID ->
   molecule = archive.get(UID)
   table = molecule.getDataTable()

   //If the column already exists we don't need to add it
   //instead we will just be overwriting the values below..
   if (!table.hasColumn("Flowrate (ul/min)"))
      table.appendColumn("Flowrate (ul/min)")

   for (int row=0;row<table.getRowCount();row++) {
      table.set("Flowrate (ul/min)", row, sliceToFlowrate.get(table.get("slice", row)))
   }

   archive.put(molecule)
}

archive.unlock()